# Cut pool cleanup for the capacity cuts in the Coluna master (PROTOTYPE).
#
# Runs through Coluna's BeforeCutGenAlgo: after each column generation and before each
# cut round (Coluna/src/Algorithm/conquer.jl, run_colcutgen!), at the root and in
# strong-branching phase 2.
#
# For every user cut (duty MasterUserCutConstr) it keeps an age = number of consecutive
# calls in which the cut was not binding (slack > CutSlackTol in the master LP primal
# solution). With CutCleanupDeactivate = true, cuts with age >= CutMinAge are
# deactivated (MathProg.deactivate!: removed from the CPLEX model at the next sync).
# Coluna records the active cuts per node (MasterCutsUnit), so a cut deactivated at a
# node stays off only in that node's subtree. A removed cut that becomes violated again
# is returned by CVRPSEP and added again as a new constraint.
#
# Prototype checks (printed when CutCleanupPrint = true):
#  - slack from the primal solution vs dual value != 0 (is the dual of the last LP
#    available here?);
#  - number of rows in the CPLEX master (does deactivate! really remove the row?).

const CutCleanupDeactivate = true    # false = only print statistics
const CutCleanupPrint      = false
const CutMinAge            = 5       # calls without being binding before removal
const CutSlackTol          = 1e-6
const CutDualTol           = 1e-9

const cutAge = Dict{Coluna.MathProg.ConstrId, Int}()   # age of each cut, global

# markRootDone = true: also sets rootDone[] (replaces RootDoneMarker in strong branching).
struct CutCleanup <: Coluna.Algorithm.AbstractOptimizationAlgorithm
    markRootDone::Bool
end

function Coluna.Algorithm.get_units_usage(::CutCleanup, reform::Coluna.MathProg.Reformulation)
    master = Coluna.MathProg.getmaster(reform)
    return [(master, Coluna.Algorithm.MasterCutsUnit, Coluna.ColunaBase.READ_AND_WRITE),
            (master, ArcFixUnit,                      Coluna.ColunaBase.READ_AND_WRITE)]
end

# Number of linear rows in the CPLEX model of the master (for the checks above).
function cplexNumRows(master)
    inner = Coluna.MathProg.getinner(Coluna.MathProg.getoptimizer(master, 1))
    F = MOI.ScalarAffineFunction{Float64}
    n = 0
    for S in (MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.EqualTo{Float64})
        n += MOI.get(inner, MOI.NumberOfConstraints{F, S}())
    end
    return n
end

function Coluna.Algorithm.run!(algo::CutCleanup, env::Coluna.Env,
                               reform::Coluna.MathProg.Reformulation,
                               conquer_output::Coluna.Algorithm.OptimizationState)
    algo.markRootDone && (rootDone[] = true)

    MP     = Coluna.MathProg
    master = MP.getmaster(reform)
    primal = Coluna.Algorithm.get_best_lp_primal_sol(conquer_output)
    isnothing(primal) && return false

    nodeArcFixing!(conquer_output)   # per-node arc fixing (NodeArcFixing.jl)

    xVal = Dict{MP.VarId, Float64}()
    for (varid, val) in primal
        xVal[varid] = val
    end

    dual = Coluna.Algorithm.get_best_lp_dual_sol(conquer_output)
    yVal = Dict{MP.ConstrId, Float64}()
    if !isnothing(dual)
        for (constrid, val) in dual
            yVal[constrid] = val
        end
    end

    # Column cleanup (ColumnCleanup.jl): needs the duals of this LP.
    if !isnothing(dual)
        zLP = Coluna.ColunaBase.getvalue(Coluna.Algorithm.get_lp_primal_bound(conquer_output))
        columnCleanup!(master, xVal, yVal, zLP)
    end

    coefmatrix = MP.getcoefmatrix(master)
    nActive = 0; nBinding = 0; nDualNz = 0; nMismatch = 0; nRemoved = 0; nInactive = 0

    for (constrid, constr) in MP.getconstrs(master)
        MP.getduty(constrid) <= MP.MasterUserCutConstr || continue
        if !MP.iscuractive(master, constr)
            nInactive += 1
            continue
        end
        nActive += 1

        lhs = 0.0
        for (varid, coef) in @view coefmatrix[constrid, :]
            lhs += coef * get(xVal, varid, 0.0)
        end
        slack   = abs(lhs - MP.getcurrhs(master, constr))
        binding = slack <= CutSlackTol
        dualNz  = abs(get(yVal, constrid, 0.0)) > CutDualTol
        nBinding += binding
        nDualNz  += dualNz
        nMismatch += (dualNz && !binding)       # nonzero dual must imply binding

        # Rank-1 cuts (Rank1Cuts.jl) take slots of the C++ bitset (MaxActiveCuts): they
        # age when their dual is 0 and are removed after R1CutMinAge calls, so the slots
        # go to new cuts. Needs the duals of this LP (no dual solution -> no aging).
        isR1   = constr.custom_data isa R1CutData
        weak   = isR1 ? (!isnothing(dual) && !dualNz) : !binding
        minAge = isR1 ? R1CutMinAge : CutMinAge

        age = weak ? get(cutAge, constrid, 0) + 1 : 0
        cutAge[constrid] = age

        if CutCleanupDeactivate && age >= minAge
            MP.deactivate!(master, constr)
            delete!(cutAge, constrid)
            nRemoved += 1
        end
    end

    if CutCleanupPrint
        println("CUTCLEANUP active $nActive | binding $nBinding | dual!=0 $nDualNz",
                " | dual!=0 but slack $nMismatch | inactive $nInactive | removed $nRemoved",
                " | cplex rows $(cplexNumRows(master)) | dual sol ", isnothing(dual) ? "none" : "ok")
    end
    return false            # no new column generation only because of the cleanup
end
