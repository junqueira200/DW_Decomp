# Up-to-date sets of the master's variables, so Coluna stops looping over every column
# ever created (Coluna 0.8.2).
#
# Profile of an E-n76-k14 tree run (2026-11-15): node records (~10%) and reading the LP
# solution (~12%) loop over ALL variables of the master, i.e. every column ever generated
# (~10^5, most of them inactive after column cleanup and node switches).
#
# Registry per formulation (kept by hooks on Coluna's only entry points):
#  - noncol: every variable that is not a column (arc representatives, artificial vars,
#    ...), in creation order; a fixed, small list;
#  - active_cols: the ids of the columns (MasterCol) that are currently active.
# Hooks: _addvar!(manager, var) (the only place a variable enters a formulation),
# _activate! / _deactivate! (the only places the active flag changes).
#
# Users of the registry (same results as Coluna's loops; only the visiting order of the
# variables changes):
#  - MasterColumnsUnit record / restore (record = copy of active_cols; restore = the
#    difference between active_cols and the record);
#  - StaticVarConstrUnit record / restore: its variable loop only concerns static duties,
#    never columns -> noncol;
#  - get_primal_solutions and LpDuals.jl's get_dual_solutions_cplex: only active
#    explicit variables matter -> noncol + active_cols.
# ActiveSetsVerify[] = true checks every use against Coluna's full scan (slow: tests only).
#
# It uses Coluna internals: check it again if Coluna is updated. Julia warns that the
# methods are overwritten.

Core.eval(Coluna.MathProg, quote

mutable struct VarRegistry
    noncol::Vector{Tuple{VarId,Variable}}
    active_cols::Set{VarId}
end
const var_registries = IdDict{Any,VarRegistry}()
const ActiveSetsVerify = Ref(false)

var_registry(m) = get!(() -> VarRegistry(Tuple{VarId,Variable}[], Set{VarId}()), var_registries, m)
is_col_id(id) = getduty(id) <= MasterCol

# Every variable that can be active: the non-columns and the active columns.
function maybe_active_vars(form)
    reg = var_registry(form.manager)
    vars = Vector{Tuple{VarId,Variable}}(undef, length(reg.noncol) + length(reg.active_cols))
    k = 0
    for p in reg.noncol
        vars[k += 1] = p
    end
    for varid in reg.active_cols
        vars[k += 1] = (varid, getvar(form, varid))
    end
    return vars
end

# Coluna's _addvar!(manager, var) (manager.jl), plus registration.
function _addvar!(m::FormulationManager, var::Variable)
    if haskey(m.vars, var.id)
        error(string(
            "Variable of id ", var.id, " exists. Its name is ", m.vars[var.id].name,
            " and you want to add a variable named ", var.name, "."
        ))
    end
    m.vars[var.id] = var
    reg = var_registry(m)
    if is_col_id(var.id)
        var.curdata.is_active && push!(reg.active_cols, var.id)
    else
        push!(reg.noncol, (var.id, var))
    end
    return
end

# Coluna's _activate! / _deactivate! (varconstr.jl), plus the active column set.
function _activate!(form::Formulation, varconstr::AbstractVarConstr)
    if isexplicit(form, varconstr) && !iscuractive(form, varconstr)
        add!(form.buffer, getid(varconstr))
    end
    varconstr.curdata.is_active = true
    if varconstr isa Variable && is_col_id(getid(varconstr))
        push!(var_registry(form.manager).active_cols, getid(varconstr))
    end
    return
end

function _deactivate!(form::Formulation, varconstr::AbstractVarConstr)
    if isexplicit(form, varconstr) && iscuractive(form, varconstr)
        remove!(form.buffer, getid(varconstr))
    end
    varconstr.curdata.is_active = false
    if varconstr isa Variable && is_col_id(getid(varconstr))
        delete!(var_registry(form.manager).active_cols, getid(varconstr))
    end
    return
end

function verify_active_sets(form)
    reg = var_registry(form.manager)
    cols = Set{VarId}()
    nnoncol = 0
    for (varid, var) in getvars(form)
        if is_col_id(varid)
            iscuractive(form, var) && push!(cols, varid)
        else
            nnoncol += 1
        end
    end
    cols == reg.active_cols || error("ActiveSets: active column set differs from the scan")
    nnoncol == length(reg.noncol) || error("ActiveSets: non-column list differs from the scan")
    return
end

# Coluna's get_primal_solutions (MOIinterface.jl), looping over the registry.
function get_primal_solutions(form::F, optimizer::MoiOptimizer) where {F <: Formulation}
    ActiveSetsVerify[] && verify_active_sets(form)
    inner = getinner(optimizer)
    nb_primal_sols = MOI.get(inner, MOI.ResultCount())
    solutions = PrimalSolution{F}[]
    for res_idx in 1:nb_primal_sols
        if MOI.get(inner, MOI.PrimalStatus(res_idx)) != MOI.FEASIBLE_POINT
            continue
        end

        solcost = getobjconst(form)
        solvars = VarId[]
        solvals = Float64[]

        for (id, var) in maybe_active_vars(form)
            iscuractive(form, var) && isexplicit(form, var) || continue
            moirec = getmoirecord(var)
            moi_index = getmoiindex(moirec)
            val = MOI.get(inner, MOI.VariablePrimal(res_idx), moi_index)
            solcost += val * getcurcost(form, id)
            val = round(val, digits = Coluna.TOL_DIGITS)
            if abs(val) > Coluna.TOL
                push!(solvars, id)
                push!(solvals, val)
            end
        end
        fixed_obj = 0.0
        for (var_id, fixed_val) in getpartialsol(form)
            push!(solvars, var_id)
            push!(solvals, fixed_val)
            fixed_obj += getcurcost(form, var_id) * fixed_val
        end
        solcost += fixed_obj
        push!(solutions, PrimalSolution(form, solvars, solvals, solcost, FEASIBLE_SOL))
    end
    return solutions
end

end)  # Core.eval(Coluna.MathProg)

Core.eval(Coluna.Algorithm, quote

# Coluna's MasterColumnsUnit record (formstorages.jl), from the active column set.
function ClB.record(::Type{MasterColumnsRecord}, id::Int, form::Formulation, unit::MasterColumnsUnit)
    MathProg.ActiveSetsVerify[] && MathProg.verify_active_sets(form)
    record = MasterColumnsRecord()
    for varid in MathProg.var_registry(form.manager).active_cols
        isexplicit(form, getvar(form, varid)) && push!(record, varid)
    end
    return record
end

# Coluna's MasterColumnsUnit restore: same end state (an explicit column is active iff it
# is in the record), reached through the difference with the current active set.
function ClB.restore_from_record!(form::Formulation, ::MasterColumnsUnit, state::MasterColumnsRecord)
    reg = MathProg.var_registry(form.manager)
    to_deactivate = VarId[]
    for varid in reg.active_cols
        if !(varid in state.active_cols) && isexplicit(form, getvar(form, varid))
            push!(to_deactivate, varid)
        end
    end
    for varid in to_deactivate
        deactivate!(form, getvar(form, varid))
    end
    for varid in state.active_cols
        var = getvar(form, varid)
        if isexplicit(form, var) && !iscuractive(form, var)
            activate!(form, var)
        end
    end
    MathProg.ActiveSetsVerify[] && MathProg.verify_active_sets(form)
    return
end

# Coluna's StaticVarConstrUnit record / restore: same code, but the variable loop goes
# over the non-column variables only (columns never have a static duty).
function ClB.record(::Type{StaticVarConstrRecord}, id::Int, form::Formulation, unit::StaticVarConstrUnit)
    record = StaticVarConstrRecord(Dict{ConstrId,ConstrState}(), Dict{VarId,StaticVarState}())
    for (id, constr) in getconstrs(form)
        if isaStaticDuty(getduty(id)) && iscuractive(form, constr) && isexplicit(form, constr)
            record.constrs[id] = ConstrState(getcurrhs(form, constr))
        end
    end
    for (id, var) in MathProg.var_registry(form.manager).noncol
        if isaStaticDuty(getduty(id)) && iscuractive(form, var) && (isexplicit(form, var) || isaMasterRepDuty(getduty(id)))
            record.vars[id] = StaticVarState(getcurcost(form, var), getcurlb(form, var), getcurub(form, var))
        end
    end
    return record
end

function ClB.restore_from_record!(form::Formulation, ::StaticVarConstrUnit, record::StaticVarConstrRecord)
    for (id, constr) in getconstrs(form)
        if isaStaticDuty(getduty(id)) && isexplicit(form, constr)
            if haskey(record.constrs, id)
                if !iscuractive(form, constr)
                    activate!(form, constr)
                end
                apply_state!(form, constr, record.constrs[id])
            else
                if iscuractive(form, constr)
                    deactivate!(form, constr)
                end
            end
        end
    end
    for (id, var) in MathProg.var_registry(form.manager).noncol
        if isaStaticDuty(getduty(id)) && (isexplicit(form, var) || isaMasterRepDuty(getduty(id)))
            if haskey(record.vars, id)
                if !iscuractive(form, var)
                    activate!(form, var)
                end
                apply_state!(form, var, record.vars[id])
            else
                if iscuractive(form, var)
                    deactivate!(form, var)
                end
            end
        end
    end
    return
end

end)  # Core.eval(Coluna.Algorithm)

# Use these from model.jl.
activeSetsReset!() = (empty!(Coluna.MathProg.var_registries); nothing)
activeSetsVerify!(on::Bool) = (Coluna.MathProg.ActiveSetsVerify[] = on; nothing)
