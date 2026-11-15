# Batch deletion of variables and constraints when Coluna syncs the master with CPLEX
# (Coluna 0.8.2 + CPLEX.jl 1.2.0).
#
# Profile of an E-n76-k14 tree run (2026-11-15): ~7% of the time in sync_solver!, almost
# all in deleting the columns (and cuts) removed since the last solve, ONE AT A TIME.
# In CPLEX.jl a single MOI.delete(model, ::VariableIndex) calls CPXdelcols for one column
# and then shifts the column index of every other variable: O(n) each, O(m*n) for m
# deletions. Column cleanup and node switching delete many columns at once.
#
# This file redefines remove_from_optimizer!(form, optimizer, ids::Set) inside
# Coluna.MathProg (called by sync_solver! for the removed constraints, then the removed
# variables). When the inner optimizer is a CPLEX.Optimizer:
#  - variables: one MOI.delete(model, Vector{VariableIndex}) (CPLEX.jl's batch method:
#    one CPXdelcols per contiguous range, one re-indexing pass). The variable's bound
#    constraints are not deleted one by one first (Coluna's code does that): in CPLEX.jl
#    the bounds belong to the variable's info and disappear with it. An integrality
#    constraint (only during a restricted master MIP) is still deleted first, as Coluna
#    does.
#  - linear constraints: CPXdelrows per contiguous range of rows and one re-indexing
#    pass of CPLEX.jl's affine_constraint_info (CPLEX.jl has no batch method for them;
#    this mirrors its single delete, MOI_wrapper.jl).
#  - Coluna's own records are reset exactly as in its single-element code.
# Otherwise (another solver): Coluna's original loop.
#
# It uses Coluna and CPLEX.jl internals (_info, _intervalize, affine_constraint_info):
# check it again if either package is updated. Julia warns that the method is overwritten.

Core.eval(Coluna.MathProg, quote

const cplex_sync_stats = zeros(Int, 4)   # batch var calls, vars deleted, batch constr calls, constrs deleted

function remove_from_optimizer!(form::Formulation, optimizer::MoiOptimizer, ids::Set{I}) where {I<:Id}
    inner = getinner(optimizer)
    if !(inner isa $(CPLEX).Optimizer) || isempty(ids)
        for id in ids                                   # Coluna's original loop
            elem = getelem(form, id)
            if elem !== nothing
                remove_from_optimizer!(form, optimizer, getelem(form, id))
            else
                definitive_deletion_from_optimizer!(form, optimizer, id)
            end
        end
        return
    end
    elems = Any[]
    for id in ids
        elem = getelem(form, id)
        if elem === nothing
            elem = id isa VarId ? form.buffer.var_buffer.definitive_deletion[id] :
                                  form.buffer.constr_buffer.definitive_deletion[id]
        end
        push!(elems, elem)
    end
    batch_remove_from_cplex!(inner, elems)
    return
end

function batch_remove_from_cplex!(inner, elems)
    vars = MOI.VariableIndex[]
    rows = MOI.ConstraintIndex[]
    for elem in elems
        moirecord = getmoirecord(elem)
        if elem isa Variable
            @assert getmoiindex(moirecord).value != -1
            if getkind(moirecord).value != -1             # integrality: delete it first
                MOI.delete(inner, getkind(moirecord))
            end
            push!(vars, getmoiindex(moirecord))
            setlowerbound!(moirecord, MoiVarLowerBound())
            setupperbound!(moirecord, MoiVarUpperBound())
            setkind!(moirecord, MoiVarKind())
            setmoiindex!(moirecord, MoiVarIndex())
        else
            @assert getmoiindex(moirecord).value != -1
            c = getmoiindex(moirecord)
            if c isa MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}}
                push!(rows, c)
            else
                MOI.delete(inner, c)                     # other kinds: one at a time
            end
            setmoiindex!(moirecord, MoiConstrIndex())
        end
    end
    if !isempty(rows)
        batch_delete_affine_rows!(inner, rows)
        cplex_sync_stats[3] += 1
        cplex_sync_stats[4] += length(rows)
    end
    if !isempty(vars)
        MOI.delete(inner, vars)                           # CPLEX.jl batch method
        cplex_sync_stats[1] += 1
        cplex_sync_stats[2] += length(vars)
    end
    return
end

# Same effect as calling CPLEX.jl's MOI.delete on each affine constraint of `cs`.
function batch_delete_affine_rows!(model, cs)
    CP = $(CPLEX)
    sorted_rows = sort!([CP._info(model, c).row for c in cs])
    starts, ends = CP._intervalize(sorted_rows)
    for ri in reverse(1:length(starts))
        ret = CP.CPXdelrows(model.env, model.lp, Cint(starts[ri] - 1), Cint(ends[ri] - 1))
        CP._check_ret(model, ret)
    end
    for c in cs
        delete!(model.affine_constraint_info, c.value)
    end
    for (key, info) in model.affine_constraint_info
        info.row -= searchsortedlast(sorted_rows, info.row)
    end
    model.name_to_constraint_index = nothing
    return
end

end)  # Core.eval(Coluna.MathProg)

# Use these from model.jl.
cplexSyncReset!() = fill!(Coluna.MathProg.cplex_sync_stats, 0)
function cplexSyncPrint()
    s = Coluna.MathProg.cplex_sync_stats
    println("CPLEX batch deletes: ", s[2], " variables in ", s[1], " calls, ",
            s[4], " constraints in ", s[3], " calls")
end
