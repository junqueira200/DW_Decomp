# Column cleanup in the master (Coluna 0.8.2 never removes columns: its cleanup_threshold
# / cleanup_ratio parameters are unused).
#
# Runs in CutCleanup.run! (after every column generation), at the root and in every
# strong-branching child, with the node's LP solution, duals and LP value.
#  (a) exact rule: a column with reduced cost > UB - zLP + (K-1)*FixRcTol cannot be in a
#      solution better than UB in this subtree (same argument as the arc fixing) ->
#      deactivated;
#  (b) size rule: above ColMaxActive active columns, the non-basic columns with the
#      largest positive reduced costs are deactivated down to ColKeepRatio*ColMaxActive.
# Per node: Coluna's MasterColumnsUnit records the active columns of every node and
# deactivates, when a node is processed, the columns not in its record, so a column
# removed here stays off only in this subtree. A deactivated column that pricing finds
# again is reactivated by Coluna (insert_columns!, colgen/default.jl).

const ColCleanupOn  = true
# Exact rule removes a column only if rc > ColFixFactor * gap. 1.0 = every provably useless
# column; > 1 keeps the borderline ones, which pricing would often generate again.
const ColFixFactor  = 2.0
const ColMaxActive  = 3000    # size rule: above this many active columns ...
const ColKeepRatio  = 0.75    # ... deactivate down to this fraction
const nColFixed     = Ref(0)  # removed by the exact rule (rc > gap)
const nColCleaned   = Ref(0)  # removed by the size rule
const maxColActive  = Ref(0)  # largest number of active columns seen
const nColNegRc     = Ref(0)  # sanity check: active columns with rc < -1e-4 (should be ~0)

function columnCleanup!(master, xVal, yVal, zLP)
    ColCleanupOn || return
    MP = Coluna.MathProg
    coefmatrix = MP.getcoefmatrix(master)
    UB  = ubArcFixing[]
    gap = (isfinite(UB) && isfinite(zLP) && zLP < UB) ?
          UB - zLP + (inst.numberOfVehicles - 1) * FixRcTol + 1e-6 : Inf

    cand = Tuple{Float64, MP.VarId}[]           # non-basic columns: (rc, id)
    nActive = 0
    for (varid, var) in MP.getvars(master)
        MP.getduty(varid) <= MP.MasterCol || continue
        (MP.iscuractive(master, var) && MP.isexplicit(master, var)) || continue
        nActive += 1
        get(xVal, varid, 0.0) > 1e-9 && continue            # in the LP solution: keep
        rc = MP.getcurcost(master, varid)
        for (constrid, coef) in @view coefmatrix[:, varid]
            rc -= coef * get(yVal, constrid, 0.0)
        end
        rc < -1e-4 && (nColNegRc[] += 1)
        if rc > ColFixFactor * gap                         # (a) exact: never useful here
            MP.deactivate!(master, var)
            nColFixed[] += 1
            nActive -= 1
        else
            push!(cand, (rc, varid))
        end
    end
    maxColActive[] = max(maxColActive[], nActive)
    if nActive > ColMaxActive                              # (b) size rule
        sort!(cand; by = first, rev = true)                # largest rc first
        nRemove = min(length(cand), nActive - round(Int, ColKeepRatio * ColMaxActive))
        for k in 1:nRemove
            rc, varid = cand[k]
            rc > 1e-6 || break                             # never remove rc <= 0 columns
            MP.deactivate!(master, varid)
            nColCleaned[] += 1
        end
    end
    return
end

columnCleanupReset!() = (nColFixed[] = 0; nColCleaned[] = 0; maxColActive[] = 0; nColNegRc[] = 0; nothing)
columnCleanupPrint() = println("Column cleanup: ", nColFixed[], " fixed (rc > gap), ", nColCleaned[],
    " cleaned (size rule), at most ", maxColActive[], " active columns, ", nColNegRc[],
    " checks with rc < -1e-4")
