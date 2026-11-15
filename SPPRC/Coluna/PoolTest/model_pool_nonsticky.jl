using JuMP, BlockDecomposition, Coluna, CPLEX
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using Coluna.Branching
using CVRPSEP

# CVRPSEP.jl moves every cut it returns into this manager's list (it grows with every
# call, duplicates included), but that list does not stop a known cut from being
# returned again when it is violated (checked 2026-11-14).
const capCutManager = Ref{CutManager}()

# Print the column generation iteration line only every ColGenPrintEvery iterations.
# Coluna 0.8.2 ignores ColumnGeneration(log_print_frequency = ...) and prints every
# iteration, so this redefines its printer method (Algorithm/colgen/printer.jl:189).
# It uses Coluna internals (_colgen_iter_str, elapsed_optim_time): check it again if
# Coluna is updated. Julia may warn that the method is overwritten.
const ColGenPrintEvery = 10

function Coluna.ColGen.after_colgen_iteration(
    ctx::Coluna.Algorithm.ColGenPrinterContext, phase, stage, env,
    colgen_iteration, stab, ip_primal_sol, colgen_iter_output)

    colgen_iteration % ColGenPrintEvery == 0 || return

    println(Coluna.Algorithm._colgen_iter_str(
        colgen_iteration, colgen_iter_output, ctx.phase,
        Coluna.ColGen.stage_id(stage), ctx.sp_elapsed_time, ctx.mst_elapsed_time,
        Coluna.elapsed_optim_time(env), Coluna.ColGen.get_output_str(stab)))
    return
end

# Time limit of the branch-and-bound (wall clock, from optimize!). Inf = none.
# Coluna 0.8.2's own time limit (TreeSearchAlgorithm(timelimit = ...) or JuMP's
# set_time_limit_sec) crashes: a node conquered after the limit has no LP solution, but
# the branching still runs on it (error("Error") in candidates_selection). So the limit
# is checked here instead, in the tree's stopping test, which runs before each node is
# taken: the search stops cleanly between nodes (it may exceed the limit by one node,
# including its strong branching). The global dual bound is still valid.
const TimeLimitSec = 600.0
const tStartOptimize = Ref(0.0)

# Redefines Coluna's stopping test of the branch-and-bound
# (Algorithm/treesearch/branch_and_bound.jl:155) to add the time limit. It uses Coluna
# internals (_update_global_dual_bound!, BaBSearchSpace fields): check it again if Coluna
# is updated. Julia may warn that the method is overwritten.
function Coluna.TreeSearch.stop(space::Coluna.Algorithm.BaBSearchSpace, untreated_nodes)
    Coluna.Algorithm._update_global_dual_bound!(space, space.reformulation, untreated_nodes)
    space.nb_untreated_nodes = length(untreated_nodes)
    timeUp = !isempty(untreated_nodes) && (time() - tStartOptimize[]) >= TimeLimitSec
    timeUp && println("Time limit of ", TimeLimitSec, " s reached: stopping the tree search.")
    return space.nb_nodes_treated >= space.max_num_nodes ||
           space.nb_untreated_nodes > space.open_nodes_limit ||
           timeUp
end

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

# Reduced-cost arc fixing, only at the root (the cut callback does not know the node).
# An arc is removed for the whole tree if every route through it has
#   z_LP + rc(route) > UB.
const ArcFixing      = true
const MaxNbCutRounds = 20       # also passed to ColCutGenConquer (root detection uses it)
const FixRcTol       = 1e-4     # same as the column submission threshold
# true  = fix only in the last root round (first round with no cut, or MaxNbCutRounds):
#         one fixing run, at the highest z_LP, but it does not help the root pricing.
# false = fix in every root round: removed arcs also shrink the graph of the next
#         root rounds, and the union over the rounds can remove more arcs.
global FixOnlyLastRound::Bool = false
const ubArcFixing    = Ref(Inf) # set from primalBound in createModel
const rootCutCalls   = Ref(0)   # cut callback calls at the root
const rootDone       = Ref(false)
const tFixing        = Ref(0.0)
#include("instanceCVRP.jl")


# Strong branching (false = Coluna's default: most fractional arc, no evaluation).
# Candidates are arcs x[1,(i,j)] (SingleVarBranchingRule), first chosen by fractionality.
#  - phase 1: SbNumCandLP candidates, each child evaluated by the restricted master LP
#             (no pricing: cheap, but optimistic);
#  - phase 2: the SbNumCandCG best ones, each child evaluated by column generation
#             and SbCutRounds capacity cut rounds (no heuristic); the best one
#             (ProductScore) is chosen.
# Coluna does NOT conquer the chosen children again: a node created by strong branching
# already has a conquer output, so the tree goes straight to its branching
# (TreeSearch.children, Algorithm/treesearch.jl:184). The phase-2 conquer is therefore the
# real conquer of every non-root node: only SbCutRounds cut rounds and no primal
# heuristic there, and these
# nodes are not counted by maxnumnodes (after_conquer! is not called for them).
global UseStrongBranching::Bool = true
const SbNumCandLP = 20
const SbNumCandCG = 3
const SbVerbose   = false   # true = print the score of every candidate
# Capacity cut rounds in strong-branching phase 2. Coluna does not conquer the chosen
# children again, so these are the only cut rounds at the non-root nodes (0 = cuts only
# at the root).
const SbCutRounds = 2

# Marks the end of the root. Runs in strong-branching phase 2 before each cut round, so
# the cut callback never treats a call from a child node as a root round (the root arc
# fixing would otherwise use a child's LP value).
struct RootDoneMarker <: Coluna.Algorithm.AbstractOptimizationAlgorithm end

function Coluna.Algorithm.get_units_usage(::RootDoneMarker, reform::Coluna.MathProg.Reformulation)
    return [(sp, Coluna.Algorithm.StaticVarConstrUnit, Coluna.ColunaBase.READ_ONLY)
            for (_, sp) in Coluna.MathProg.get_dw_pricing_sps(reform)]
end

function Coluna.Algorithm.run!(::RootDoneMarker, env::Coluna.Env,
                               reform::Coluna.MathProg.Reformulation,
                               conquer_output::Coluna.Algorithm.OptimizationState)
    rootDone[] = true
    return false                                   # the master is not changed
end

include(joinpath(@__DIR__, "CutCleanup.jl"))   # cut pool cleanup (CutCleanup)

# Same column generation everywhere (tree nodes and strong branching).
colgenAlgo() = Coluna.Algorithm.ColumnGeneration(
    max_nb_iterations = 1000000,
    #show_column_already_inserted_warning = true,
    smoothing_stabilization = 0.5,
    #stages_pricing_solver_ids = [2, 1],
    redcost_tol = 1E-3
)

function strongBranchingAlgo()
    return Coluna.Algorithm.StrongBranching(
        phases = [
            Coluna.Algorithm.BranchingPhase(
                SbNumCandLP,
                Coluna.Algorithm.RestrMasterLPConquer(),
                Coluna.Algorithm.ProductScore()),
            Coluna.Algorithm.BranchingPhase(
                SbNumCandCG,
                Coluna.Algorithm.ColCutGenConquer(
                    colgen = colgenAlgo(),
                    primal_heuristics = Coluna.Algorithm.ParameterizedHeuristic[],
                    before_cutgen_user_algorithm =
                        Coluna.Algorithm.BeforeCutGenAlgo(CutCleanup(true), "Cut cleanup + root done"),
                    max_nb_cut_rounds = SbCutRounds),
                Coluna.Algorithm.ProductScore()),
        ],
        rules = [Coluna.Branching.PrioritisedBranchingRule(
                    Coluna.Algorithm.SingleVarBranchingRule(), 1.0, 1.0)],
        selection_criterion = Coluna.Algorithm.MostFractionalCriterion(),
        verbose = SbVerbose
    )
end

function createModel(inst::CvrpData; primalBound::Float64 = Inf)


    @axis(VehicleTypes, [1])

    dividealg = UseStrongBranching ? strongBranchingAlgo() : Coluna.Algorithm.ClassicBranching()

    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
        solver = Coluna.Algorithm.TreeSearchAlgorithm(
            #opennodeslimit = 100,
            maxnumnodes = 50,
            explorestrategy = Coluna.TreeSearch.BestDualBoundStrategy(),
            dividealg = dividealg,
            conqueralg = Coluna.Algorithm.ColCutGenConquer(
                colgen = colgenAlgo(),
                # Restricted master heuristic only at the root (max_depth = 0).
                # Coluna's default runs it at every node (its `frequency` field is not
                # used in Coluna 0.8.2), solving a MIP each time.
                # Arguments: algorithm, root_priority, nonroot_priority, frequency,
                # max_depth, name. Use ParameterizedHeuristic[] to turn it off.
                primal_heuristics = [Coluna.Algorithm.ParameterizedHeuristic(
                    Coluna.Algorithm.RestrictedMasterHeuristic(),
                    1.0, 1.0, 1, 0, "Restricted Master IP")],
                before_cutgen_user_algorithm =
                    Coluna.Algorithm.BeforeCutGenAlgo(CutCleanup(false), "Cut cleanup"),
                # Cut rounds per node (Coluna's default is 3; the root needs more).
                max_nb_cut_rounds = MaxNbCutRounds
            )
        )
    ),
        "default_optimizer" => optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPXPARAM_ScreenOutput" => 0,
            # Master LPs are small and re-solved many times: thread start-up costs
            # more than parallelism gains. TEST: compare with the default (0 = auto).
            "CPXPARAM_Threads" => 1,
            # Primal simplex: adding columns keeps the previous basis primal feasible,
            # so re-solving the master starts close to the optimum.
            # TEST: 0 = automatic (default), 1 = primal, 2 = dual simplex.
            "CPXPARAM_LPMethod" => 1
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
    
    # Known upper bound (e.g. optimum + 1) so nodes can be pruned; Inf = none.
    if isfinite(primalBound)
        objectiveprimalbound!(model, primalBound)
    end

    # Arc fixing state (root only).
    ubArcFixing[]  = primalBound
    rootCutCalls[] = 0
    rootDone[]     = false
    empty!(cutAge)                 # cut ages (CutCleanup.jl)
    empty!(cutKeyToId); empty!(cutSticky); nCutReact[] = 0; nCutSkipped[] = 0
    fill!(removedArcs, false)

    # Rounded capacity cuts (CVRPSEP), separated after column generation at each node.
    capCutManager[] = CutManager()
    MOI.set(model, MOI.UserCutCallback(),
            cbdata -> capacityCutCallback(cbdata, mdl, xVar, inst))

    tStartOptimize[] = time()
    optimize!(model)

    println("Time in solveSpprc: ", tSpprc[], " s")
    println("  heuristic: ", tHeur[],  " s in ", nHeur[],  " calls")
    println("  exact:     ", tExact[], " s in ", nExact[], " calls")
    println("Time in arc fixing: ", tFixing[], " s; arcs removed: ", count(removedArcs))
    println("Cut pool: reactivated ", nCutReact[], ", skipped (already in pool) ", nCutSkipped[])

end

# Rounded capacity cuts: x(S:S) <= |S| - k(S), separated by CVRPSEP.
# CVRPSEP numbering: depot = 1, customers = 2..n  (ours: depot = 0, customers = 1..n-1).
function capacityCutCallback(cbdata, model, x, inst)
    n = inst.n

    # Arc values -> edge values x(i,j) + x(j,i), only positive edges.
    # zLP = value of the master LP (sum of c(i,j) * x(i,j)), used by the arc fixing.
    edgeTail = Int[]; edgeHead = Int[]; edgeX = Float64[]
    zLP = 0.0
    for i in 0:(n - 1), j in (i + 1):(n - 1)
        vij = callback_value(cbdata, x[1, (i, j)])
        vji = callback_value(cbdata, x[1, (j, i)])
        zLP += inst.arc_weights[(i, j)] * vij + inst.arc_weights[(j, i)] * vji
        v = vij + vji
        if v > 1e-6
            push!(edgeTail, i + 1); push!(edgeHead, j + 1); push!(edgeX, v)
        end
    end

    numCuts = 0
    if !isempty(edgeX)
        demand = [inst.demand[i] for i in 0:(n - 1)]      # demand[1] = depot (ignored)

        S, RHS = rounded_capacity_cuts!(capCutManager[], demand, inst.capacity,
                                        edgeTail, edgeHead, edgeX;
                                        integrality_tolerance = 1e-6, max_n_cuts = 5)

        for (s, rhs) in zip(S, RHS)
            nodes = s .- 1                                  # back to our numbering
            key   = sort(nodes)
            if haskey(cutKeyToId, key)          # already in the master (CutCleanup
                nCutSkipped[] += 1              # reactivates it if it is violated)
                continue
            end
            cut = @build_constraint(
                sum(x[1, (i, j)] for i in nodes, j in nodes if i != j) <= rhs)
            cutKeyToId[key] = MOI.submit(model, MOI.UserCut(cbdata), cut)   # Coluna returns the ConstrId
            numCuts += 1
        end
        printInfo && println("capacity cuts added: ", numCuts)
    end

    # Root = the cut rounds until the first round with no cut, or MaxNbCutRounds rounds
    # (Coluna runs the cut callback right after column generation converges).
    if ArcFixing && !rootDone[]
        rootCutCalls[] += 1
        lastRootRound = (numCuts == 0 || rootCutCalls[] >= MaxNbCutRounds)
        if !FixOnlyLastRound || lastRootRound
            fixArcsAtRoot(zLP, inst)
        end
        rootDone[] = lastRootRound
    end
    return
end

# Reduced-cost arc fixing with the duals of the last pricing call (vetReducedCosts).
# For a solution that uses arc a:
#   cost >= zLP + arcMin(a) + (K-1)*min(0, r*),   r* = min route reduced cost
# so arc a is removed (whole tree) if arcMin(a) > UB - zLP + (K-1)*FixRcTol.
function fixArcsAtRoot(zLP, inst)
    UB = ubArcFixing[]
    isfinite(UB) || return
    n = inst.n
    K = inst.numberOfVehicles

    gap = UB - zLP + (K - 1) * FixRcTol + 1e-6

    t = @elapsed minRouteRc = ccall((:fixArcsSpprc, LIB), Cdouble,
                                    (Ptr{Cdouble}, Cdouble, Ptr{Cdouble}),
                                    vetReducedCosts, gap, vetArcMin)
    tFixing[] += t

    # vetReducedCosts must be the converged (non-smoothed) duals of this master LP.
    if minRouteRc < -FixRcTol
        println("arc fixing skipped: min route rc = ", minRouteRc, " (duals not converged)")
        return
    end

    newly = 0
    for k in eachindex(vetArcMin)
        if !removedArcs[k] && isfinite(vetReducedCosts[k]) && vetArcMin[k] > gap
            removedArcs[k] = true
            newly += 1
        end
    end
    println("arc fixing: zLP = ", round(zLP, digits = 2), ", gap = ", round(gap, digits = 2),
            ", removed ", newly, " (total ", count(removedArcs), " of ", n*(n - 1), ") in ",
            round(t, digits = 3), " s")
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

            if i == j || removedArcs[k]      # removed by arc fixing: also skips the callback
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
    isComplete = Ref{Cint}(0)           # out: 1 = exact search finished
    exact = false
    tCall = @elapsed numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint, Ref{Cint}), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                0,
                isComplete)
    tSpprc[] += tCall
    tHeur[]  += tCall
    nHeur[]  += 1

    
    if numGeneratedRoutes == 0
        printInfo && println("--> Running exact pricing")
        exact = true

        tCall = @elapsed numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint, Ref{Cint}), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                1,
                isComplete)
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

    if exact && isComplete[] == 1
        # Complete exact search: best_rc is the minimum reduced cost.
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
    isComplete = Ref{Cint}(0)           # out: 1 = exact search finished
    numGeneratedRoutes = 
        ccall((:solveSpprc, LIB), 
                Cint, 
                (Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}, Ptr{Cdouble}, Cint, Ref{Cint}), 
                vetReducedCosts,
                vetRoutesOutput,
                sizeVetRoutesOutput,
                vetReducedCostsRoutes,
                1,
                isComplete)

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

    # Complete exact search: best_rc is a valid bound; stopped early: no bound.
    MOI.submit(model, BD.PricingDualBound(cbdata), isComplete[] == 1 ? best_rc : -Inf)

    printInfo && println("\tBest reduced cost: $best_rc; numGeneratedRoutes: $numGeneratedRoutes")
    #println("\n\n********************************\n")
    return
end