using JuMP, BlockDecomposition, Coluna, CPLEX
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using Coluna.Branching

# Prints per-call information of the pricing (number of routes, best reduced cost, ...)
global printInfo::Bool = false
# Read the subproblem variable bounds (callback_lb / callback_ub) in every pricing call.
# Coluna branches on original variables through master constraints, so the subproblem
# bounds stay 0/1 (checked with branching: no arc was ever fixed). Turn it on only if
# presolve or another branching rule can change the subproblem bounds: it costs
# 2 callbacks per arc per call.
global checkSpBounds::Bool = false
# Total time spent inside solveSpprc (C++ labeling). Reset with tSpprc[] = 0.0
const tSpprc = Ref(0.0)
# Time and number of calls of solveSpprc, split by heuristic / exact labeling
const tHeur  = Ref(0.0); const nHeur  = Ref(0)
const tExact = Ref(0.0); const nExact = Ref(0)
#include("instanceCVRP.jl")


# Phase 1: Test candidates on Restricted Master LP, keep top 3
phase_fast = Coluna.Algorithm.BranchingPhase(
    3,
    Coluna.Algorithm.ColCutGenConquer(),
    Coluna.Algorithm.ProductScore()
)

# Phase 2: Run full Column Generation on top 3 candidates, keep best 1
phase_exact = Coluna.Algorithm.BranchingPhase(
    1,
    Coluna.Algorithm.ColCutGenConquer(),
    Coluna.Algorithm.ProductScore()
)

function createModel(inst::CvrpData)

    
    @axis(VehicleTypes, [1])

    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
        solver = Coluna.Algorithm.TreeSearchAlgorithm(
            opennodeslimit = 100,
            maxnumnodes = 1,
            explorestrategy = Coluna.TreeSearch.BestDualBoundStrategy(),
            conqueralg = Coluna.Algorithm.ColCutGenConquer(
                colgen = Coluna.Algorithm.ColumnGeneration(
                    max_nb_iterations=1000000,
                    #show_column_already_inserted_warning = true,
                    smoothing_stabilization = 0.85,
                    #stages_pricing_solver_ids = [2, 1],
                    redcost_tol=1E-3
                ),
                # TEST: restricted master heuristic off, to measure its cost
                #primal_heuristics = Coluna.Algorithm.ParameterizedHeuristic[]
            )
        )
    ),
        "default_optimizer" => optimizer_with_attributes(
            CPLEX.Optimizer,             
            "CPXPARAM_ScreenOutput" => 0
        ),
    )

    global model = BlockModel(coluna)
    global x = @variable(model, x[v in VehicleTypes, a in inst.arcs], Bin)
    #println(x)
    #println(inst.arcs)
    
    @constraint(model, visit[i in inst.customers], sum(x[1, a] for a in inst.arcs if a[2] == i) == 1)
    @constraint(model, numVeicles, sum(x[1, (0, i)] for i in inst.customers) <= inst.numberOfVehicles)
    @objective(model, Min, sum(inst.arc_weights[a] * x[1, a] for a in inst.arcs))

    @dantzig_wolfe_decomposition(model, decomposition, VehicleTypes)
    subproblems = BD.getsubproblems(decomposition)
    # The closure passes model and x as arguments: inside the callback their
    # concrete types are known, so there is no dynamic dispatch on globals.
    mdl  = model
    xVar = x
    BD.specify!.(subproblems,
                 lower_multiplicity = 0,
                 upper_multiplicity = length(inst.customers),
                 solver = cbdata -> pricing_callback_heuristic(cbdata, mdl, xVar))
    
    #objectiveprimalbound!(model, 685.0)
    optimize!(model)

    println("Time in solveSpprc: ", tSpprc[], " s")
    println("  heuristic: ", tHeur[],  " s in ", nHeur[],  " calls")
    println("  exact:     ", tExact[], " s in ", nExact[], " calls")

end

function linearIndex(i::Int, j::Int, n::Int)
    return i*n+j
end

function populateReducedCosts(cbdata, model, x)

    sp_id = BD.callback_spid(cbdata, model)
    n     = inst.n

    # Arcs fixed to 1 by branching. They are applied after every reduced cost has
    # been written, so step 1 cannot overwrite them.
    fixedArcs = Tuple{Int, Int}[]

    # 1. Reduced costs and arcs fixed to 0, in a single pass.
    #    createArcs() creates every arc (i, j) with i != j, so no search in
    #    inst.arcs is needed.
    for i in 0:(n - 1)
        for j in 0:(n - 1)
            k = linearIndex(i, j, n) + 1

            if i == j
                vetReducedCosts[k] = Inf
                continue
            end

            var = x[sp_id, (i, j)]
            vetReducedCosts[k] = BD.callback_reduced_cost(cbdata, var)

            if checkSpBounds
                if BD.callback_ub(cbdata, var) < 0.5      # arc fixed to 0
                    vetReducedCosts[k] = Inf
                    printInfo && println("populateReducedCosts: arc fixed to 0: ", (i, j))
                elseif BD.callback_lb(cbdata, var) > 0.5  # arc fixed to 1
                    push!(fixedArcs, (i, j))
                end
            end
        end
    end

    # 2. Arc (i, j) fixed to 1:
    #    - a route that visits customer i must leave it through j;
    #    - a route that visits customer j must enter it from i.
    #    The depot (0) is shared by every route, so its arcs are never blocked
    #    because of another arc.
    for (i, j) in fixedArcs
        if i != 0                       # i is a customer: its only exit is j
            for jj in 0:(n - 1)
                if jj != j && jj != i
                    vetReducedCosts[linearIndex(i, jj, n) + 1] = Inf
                end
            end
        end

        if j != 0                       # j is a customer: its only entry is i
            for ii in 0:(n - 1)
                if ii != i && ii != j
                    vetReducedCosts[linearIndex(ii, j, n) + 1] = Inf
                end
            end
        end
    end

    if !isempty(fixedArcs)
        printInfo && println("populateReducedCosts: arcs fixed to 1 in the subproblem: ", fixedArcs)
    end
end

function parseRoutesFromVetRoutesOutput(actual_size)
    
    # 3. Parse routes from vetRoutesOutput (e.g. 0 1 2 0 0 5 6 7 0 ... 0 9 0)
    routes = Vector{Vector{Int}}()
    current_route = Int[]
    in_route = false

    for idx in 1:actual_size
        node = vetRoutesOutput[idx]
        if !in_route
            if node == 0
                in_route = true
                current_route = [0]
            end
        else
            push!(current_route, node)
            if node == 0
                # End of a route (at least depot -> customer -> depot)
                if length(current_route) >= 3
                    push!(routes, copy(current_route))
                end
                in_route = false
            end
        end
    end

    return routes
end

function submitRoutesToColuna(routes::Vector{Vector{Int}}, cbdata, sp_id, model, x)
    printInfo && println("**********************************\n\n")
    best_rc = Inf

    if isempty(routes)
        return 0
    end
    
    for route in routes
        # Build arcs for the route: (0, 1), (1, 2), (2, 0)
        sol_vars = [x[sp_id, (route[k], route[k + 1])] for k in 1:(length(route) - 1)]
        sol_vals = fill(1.0, length(sol_vars))

        # Calculate reduced cost of the column
        sol_rc = sum(BD.callback_reduced_cost(cbdata, v) for v in sol_vars)

        if sol_rc < best_rc
            best_rc = sol_rc
        end

        # Submit column if it has a negative reduced cost (improving column)
        if sol_rc < -1e-4
            MOI.submit(model, BD.PricingSolution(cbdata), sol_rc, sol_vars, sol_vals)

        else

            printInfo && println("Error, route is not negative. Value: ", sol_rc)

        end
        #println(sol_rc)
    end

    return best_rc
end

function pricing_callback_heuristic(cbdata, model, x)

    printInfo && println("--> Running heuristic pricing callback")
    #@show propertynames(cbdata.form)
    #exit(-1)

    dual_bound = cbdata.dual_bound
    primal_solutions = cbdata.primal_solutions

    #@info "primal_solutions: " propertynames(cbdata.primal_solutions)
    #@info "cbdata fields: " propertynames(cbdata)
    #println("--> Running pricing callback")

    sp_id = BD.callback_spid(cbdata, model)
    
    populateReducedCosts(cbdata, model, x)

    numGeneratedRoutes = 0

    # 2. Call the SPPRC solver
    sizeVetRoutesOutput = Ref{Cint}(0)
    exact = false
    tCall = @elapsed numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                0)
    tSpprc[] += tCall
    tHeur[]  += tCall
    nHeur[]  += 1

    
    if numGeneratedRoutes == 0
        printInfo && println("--> Running exact pricing")
        exact = true

        tCall = @elapsed numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                1)
        tSpprc[] += tCall
        tExact[] += tCall
        nExact[] += 1

        printInfo && println("numGeneratedRoutes: ", numGeneratedRoutes)

        if numGeneratedRoutes == 0
            MOI.submit(model, BD.PricingDualBound(cbdata), 0)
            return
        end
        

    end
    

    actual_size = sizeVetRoutesOutput[]
    #println("Number of integers returned: ", actual_size)



    # If no routes were found
    if actual_size == 0
        MOI.submit(model, BD.PricingDualBound(cbdata), -Inf)
        return
    end

    # 3. Parse routes from vetRoutesOutput (e.g. 0 1 2 0 0 5 6 7 0 ... 0 9 0)
    routes = parseRoutesFromVetRoutesOutput(actual_size)

    best_rc = submitRoutesToColuna(routes, cbdata, sp_id, model, x)

    if exact
        MOI.submit(model, BD.PricingDualBound(cbdata), best_rc)
    else
        MOI.submit(model, BD.PricingDualBound(cbdata), -Inf)
    end

    printInfo && println("\tBest reduced cost: $best_rc; numGeneratedRoutes: $numGeneratedRoutes")
    #println("\n\n********************************\n")
    return
end

function pricing_callback_exact(cbdata, model, x)
    printInfo && println("--> Running exact pricing callback")

    populateReducedCosts(cbdata, model, x)
    sp_id = BD.callback_spid(cbdata, model)

    # 2. Call the SPPRC solver
    sizeVetRoutesOutput = Ref{Cint}(0)
    numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                1)

    actual_size = sizeVetRoutesOutput[]
    #println("Number of integers returned: ", actual_size)


    # If no routes were found
    if actual_size == 0
        MOI.submit(model, BD.PricingDualBound(cbdata), 0)
        return
    end

    # 3. Parse routes from vetRoutesOutput (e.g. 0 1 2 0 0 5 6 7 0 ... 0 9 0)
    routes = parseRoutesFromVetRoutesOutput(actual_size)

    best_rc = submitRoutesToColuna(routes, cbdata, sp_id, model, x)

    MOI.submit(model, BD.PricingDualBound(cbdata), best_rc)

    printInfo && println("\tBest reduced cost: $best_rc; numGeneratedRoutes: $numGeneratedRoutes")
    #println("\n\n********************************\n")
    return
end