using JuMP, BlockDecomposition, Coluna, CPLEX
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using Coluna.Branching
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
            opennodeslimit = 1000000,
            explorestrategy = Coluna.TreeSearch.BestDualBoundStrategy(),
            conqueralg = Coluna.Algorithm.ColCutGenConquer(
                colgen = Coluna.Algorithm.ColumnGeneration(
                    max_nb_iterations=1000000,
                    #show_column_already_inserted_warning = true,
                    smoothing_stabilization = 0.85,
                    stages_pricing_solver_ids = [2, 1],
                    redcost_tol=1E-3
                )
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
    BD.specify!.(subproblems,
                 lower_multiplicity = 1,
                 upper_multiplicity = inst.numberOfVehicles,
                 solver = [pricing_callback_heuristic, pricing_callback_exact])
    
    objectiveprimalbound!(model, 685.0)
    optimize!(model)

end

function linearIndex(i::Int, j::Int, n::Int)
    return i*n+j
end

function populateReducedCosts(cbdata)

    sp_id = BD.callback_spid(cbdata, model)

    # 1. Populate reduced costs for all arcs
    for i in 0:(inst.n - 1)
        for j in 0:(inst.n - 1)
            k = linearIndex(i, j, inst.n) + 1
            rc = Inf
            if i != j && (i, j) in inst.arcs
                rc = BD.callback_reduced_cost(cbdata, x[sp_id, (i, j)])
            end
            vetReducedCosts[k] = rc
        end
    end

    for i in 0:(inst.n - 1)
        for j in 0:(inst.n - 1)

            if i == 0 || !((i, j) in inst.arcs)
                continue
            end

            lb = BD.callback_lb(cbdata, x[sp_id, (i, j)])
            up = BD.callback_ub(cbdata, x[sp_id, (i, j)])

            if up == 0
                k = linearIndex(i, j, inst.n) + 1
                vetReducedCosts[k] = Inf
            end

            if lb == 1
                for ii in 0:(inst.n -1)
                    if i != ii
                        k = linearIndex(ii, j, inst.n) + 1
                        vetReducedCosts[k] = Inf
                    end
                end
            end
        end
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

function submitRoutesToColuna(routes::Vector{Vector{Int}}, cbdata, sp_id)
    println("**********************************\n\n")
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

            println("Error, route is not negative. Value: ", sol_rc)

        end
        #println(sol_rc)
    end

    return best_rc
end

function pricing_callback_heuristic(cbdata)

    println("--> Running heuristic pricing callback")
    #@show propertynames(cbdata.form)
    #exit(-1)

    dual_bound = cbdata.dual_bound
    primal_solutions = cbdata.primal_solutions

    #@info "primal_solutions: " propertynames(cbdata.primal_solutions)
    #@info "cbdata fields: " propertynames(cbdata)
    #println("--> Running pricing callback")

    sp_id = BD.callback_spid(cbdata, model)
    
    populateReducedCosts(cbdata)

    numGeneratedRoutes = 0

    # 2. Call the SPPRC solver
    sizeVetRoutesOutput = Ref{Cint}(0)
    exact = false
    numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                0)

    #=
    if numGeneratedRoutes == 0
        println("--> Running exact pricing")
        exact = true

        numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                1)

        println("numGeneratedRoutes: ", numGeneratedRoutes)

        if numGeneratedRoutes == 0
            MOI.submit(model, BD.PricingDualBound(cbdata), 0)
            return
        end
        

    end
    =#

    actual_size = sizeVetRoutesOutput[]
    #println("Number of integers returned: ", actual_size)



    # If no routes were found
    if actual_size == 0
        MOI.submit(model, BD.PricingDualBound(cbdata), -Inf)
        return
    end

    # 3. Parse routes from vetRoutesOutput (e.g. 0 1 2 0 0 5 6 7 0 ... 0 9 0)
    routes = parseRoutesFromVetRoutesOutput(actual_size)

    best_rc = submitRoutesToColuna(routes, cbdata, sp_id)

    #if exact
    #    MOI.submit(model, BD.PricingDualBound(cbdata), best_rc)
    #else
    MOI.submit(model, BD.PricingDualBound(cbdata), -Inf)
    #end

    println("\tBest reduced cost: $best_rc; numGeneratedRoutes: $numGeneratedRoutes")
    #println("\n\n********************************\n")
    return
end

function pricing_callback_exact(cbdata)
    println("--> Running exact pricing callback")

    populateReducedCosts(cbdata)
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

    best_rc = submitRoutesToColuna(routes, cbdata, sp_id)

    MOI.submit(model, BD.PricingDualBound(cbdata), best_rc)

    println("\tBest reduced cost: $best_rc; numGeneratedRoutes: $numGeneratedRoutes")
    #println("\n\n********************************\n")
    return
end