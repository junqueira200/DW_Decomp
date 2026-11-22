# Faster computation of a new cut's coefficients on the columns (Coluna 0.8.2).
#
# Profile of an E-n76-k7 tree run (2026-11-15): 11% of the time in Coluna's
# _setrobustmembers!(form, constr, members) (MathProg/formulation.jl), called when a cut
# (or a branching constraint) is added to the master. For each arc variable of the cut it
# visits every column of the pricing solution pool that uses that arc and does
#     coef_matrix[cut, column] += column_coeff * arc_coeff
# i.e. one read and one write in the sparse matrix (both storage directions) per
# (arc, column) pair; a column with several arcs in the cut is read and written several
# times.
#
# This file redefines it inside Coluna.MathProg: the column coefficients are first summed
# in a dictionary (same summation order per column: existing value, then the arcs in the
# order of `members`), then written once per column. Same matrix at the end. Every column
# of the pool still gets its coefficient (an inactive column may be reactivated later).
#
# It uses Coluna internals: check it again if Coluna is updated. Julia warns that the
# method is overwritten.

Core.eval(Coluna.MathProg, quote

const cut_insert_stats = zeros(Int, 3)   # calls, (arc, column) pairs, distinct columns written

function _setrobustmembers!(form::Formulation, constr::Constraint, members::VarMembership)
    coef_matrix = getcoefmatrix(form)
    constrid = getid(constr)
    acc = Dict{VarId,Float64}()
    npairs = 0
    for (varid, var_coeff) in members
        # Coefficient of the variable itself (as in Coluna's code).
        coef_matrix[constrid, varid] = var_coeff

        if getduty(varid) <= MasterRepPricingVar || getduty(varid) <= MasterRepPricingSetupVar
            for (_, spform) in get_dw_pricing_sps(form.parent_formulation)
                for (col_id, col_coeff) in @view get_primal_sol_pool(spform).solutions[:, varid]
                    # First visit of the column: start from its current coefficient,
                    # exactly as Coluna's first += would.
                    cur = get(acc, col_id, nothing)
                    if cur === nothing
                        cur = coef_matrix[constrid, col_id]
                    end
                    acc[col_id] = cur + col_coeff * var_coeff
                    npairs += 1
                end
            end
        end
    end
    for (col_id, val) in acc
        coef_matrix[constrid, col_id] = val
    end
    cut_insert_stats[1] += 1
    cut_insert_stats[2] += npairs
    cut_insert_stats[3] += length(acc)
    return
end

end)  # Core.eval(Coluna.MathProg)

# Use these from model.jl.
cutInsertReset!() = fill!(Coluna.MathProg.cut_insert_stats, 0)
function cutInsertPrint()
    s = Coluna.MathProg.cut_insert_stats
    println("Cut insertion: ", s[1], " constraints, ", s[2], " (arc, column) pairs, ",
            s[3], " column coefficients written")
end
