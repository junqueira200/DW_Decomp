# Faster reading of the master LP duals for Coluna 0.8.2 + CPLEX.jl 1.2.0.
#
# Profile of an E-n101 tree run (2026-11-15): Coluna's get_dual_solutions
# (Coluna/src/MathProg/MOIinterface.jl) took about as long as CPLEX solving the column
# generation master LPs. Per LP solve it
#   - calls MOI.get(ConstraintDual) once per constraint: one CPXgetpi call per row;
#   - loops over EVERY variable of the master (all columns ever generated, active or
#     not), calling MOI.is_valid up to 4 times per variable and MOI.get(ConstraintDual)
#     once per valid bound: one CPXgetdj call per column.
# Split of its time: is_valid ~38%, Coluna's loop ~32%, single CPLEX calls ~24%.
#
# This file redefines get_dual_solutions inside Coluna.MathProg (same results):
#   - when the inner optimizer is a CPLEX.Optimizer with one result and no dual
#     certificate: all row duals with one CPXgetpi and all reduced costs with one
#     CPXgetdj, then the same formulas and sign rules as CPLEX.jl's ConstraintDual
#     getters (MOI_wrapper.jl: _dual_multiplier for rows, the min/max heuristic for
#     the LessThan / GreaterThan variable bounds);
#   - inactive or non-explicit variables and inactive constraints are skipped before
#     MOI.is_valid (they are not in the CPLEX model: Coluna removes them when it syncs
#     before each solve, so is_valid would be false anyway);
#   - otherwise: Coluna's original code (copied below, unchanged).
#
# It uses Coluna and CPLEX.jl internals (CPLEX._info, CPLEX.column, has_dual_certificate):
# check it again if either package is updated. Julia warns that the method is overwritten.

Core.eval(Coluna.MathProg, quote

const lp_duals_stats = zeros(Int, 2)   # fast path calls, original path calls

# Row dual of a master constraint, from the bulk CPXgetpi vector.
lp_row_dual(inner, c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any}, pi, mult) =
    mult * pi[$(CPLEX)._info(inner, c).row]
lp_row_dual(inner, c, pi, mult) = MOI.get(inner, MOI.ConstraintDual(1), c)

# Variable bound duals from the bulk CPXgetdj vector, CPLEX.jl's rules.
function lp_bound_dual(inner, c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}, dj, is_min_moi, mult)
    p = dj[$(CPLEX).column(inner, c)]
    is_min_moi && p > 0 && return p
    !is_min_moi && p < 0 && return -p
    return 0.0
end
function lp_bound_dual(inner, c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}, dj, is_min_moi, mult)
    p = dj[$(CPLEX).column(inner, c)]
    is_min_moi && p < 0 && return p
    !is_min_moi && p > 0 && return -p
    return 0.0
end
lp_bound_dual(inner, c::MOI.ConstraintIndex{MOI.VariableIndex,<:Union{MOI.Interval{Float64},MOI.EqualTo{Float64}}}, dj, is_min_moi, mult) =
    mult * dj[$(CPLEX).column(inner, c)]
lp_bound_dual(inner, c, dj, is_min_moi, mult) = MOI.get(inner, MOI.ConstraintDual(1), c)

function get_dual_solutions(form::F, optimizer::MoiOptimizer) where {F <: Formulation}
    inner = getinner(optimizer)
    if inner isa $(CPLEX).Optimizer && !inner.has_dual_certificate &&
       MOI.get(inner, MOI.ResultCount()) == 1
        lp_duals_stats[1] += 1
        return get_dual_solutions_cplex(form, inner)
    end
    lp_duals_stats[2] += 1
    return get_dual_solutions_original(form, optimizer)
end

function get_dual_solutions_cplex(form::F, inner) where {F <: Formulation}
    solutions = DualSolution{F}[]
    sense = getobjsense(form) == MinSense ? 1.0 : -1.0
    MOI.get(inner, MOI.DualStatus(1)) != MOI.FEASIBLE_POINT && return solutions

    env, lp = inner.env, inner.lp
    nrows = $(CPLEX).CPXgetnumrows(env, lp)
    ncols = $(CPLEX).CPXgetnumcols(env, lp)
    pi = Vector{Cdouble}(undef, nrows)
    dj = Vector{Cdouble}(undef, ncols)
    nrows > 0 && $(CPLEX)._check_ret(inner, $(CPLEX).CPXgetpi(env, lp, pi, 0, nrows - 1))
    ncols > 0 && $(CPLEX)._check_ret(inner, $(CPLEX).CPXgetdj(env, lp, dj, 0, ncols - 1))
    is_min_moi = MOI.get(inner, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    mult = is_min_moi ? 1.0 : -1.0       # CPLEX.jl's _dual_multiplier

    solcost = getobjconst(form)

    # Dual values of the constraints
    solconstrs = ConstrId[]
    solvals = Float64[]
    for (id, constr) in getconstrs(form)
        iscuractive(form, constr) || continue         # not in CPLEX
        moi_index = getmoiindex(getmoirecord(constr))
        MOI.is_valid(inner, moi_index) || continue
        val = lp_row_dual(inner, moi_index, pi, mult)
        solcost += val * getcurrhs(form, id)
        val = round(val, digits = Coluna.TOL_DIGITS)
        if abs(val) > Coluna.TOL
            push!(solconstrs, id)
            push!(solvals, sense * val)
        end
    end

    # Dual values & active bounds of the variables
    varids = VarId[]
    varvals = Float64[]
    activebounds = ActiveBound[]
    for (varid, var) in maybe_active_vars(form)                       # ActiveSets.jl
        (iscuractive(form, var) && isexplicit(form, var)) || continue   # not in CPLEX
        moirec = getmoirecord(var)
        MOI.is_valid(inner, getmoiindex(moirec)) || continue
        moi_lower_bound_index = getlowerbound(moirec)
        if MOI.is_valid(inner, moi_lower_bound_index)
            val = lp_bound_dual(inner, moi_lower_bound_index, dj, is_min_moi, mult)
            if abs(val) > Coluna.TOL
                solcost += val * getcurlb(form, varid)
                push!(varids, varid)
                push!(varvals, sense * val)
                push!(activebounds, LOWER)
            end
        end
        moi_upper_bound_index = getupperbound(moirec)
        if MOI.is_valid(inner, moi_upper_bound_index)
            val = lp_bound_dual(inner, moi_upper_bound_index, dj, is_min_moi, mult)
            if abs(val) > Coluna.TOL
                solcost += val * getcurub(form, varid)
                push!(varids, varid)
                push!(varvals, sense * val)
                push!(activebounds, UPPER)
            end
        end
    end
    fixed_obj = 0.0
    for (var_id, fixed_val) in getpartialsol(form)
        cost = getcurcost(form, var_id)
        if abs(cost) > Coluna.TOL
            fixed_obj += cost * fixed_val
        end
    end
    solcost += fixed_obj
    push!(solutions, DualSolution(
        form, solconstrs, solvals, varids, varvals, activebounds, sense*solcost,
        FEASIBLE_SOL
    ))
    return solutions
end

# Coluna's original get_dual_solutions (MOIinterface.jl:263), unchanged.
function get_dual_solutions_original(form::F, optimizer::MoiOptimizer) where {F <: Formulation}
    inner = getinner(optimizer)
    nb_dual_sols = MOI.get(inner, MOI.ResultCount())
    solutions = DualSolution{F}[]
    sense = getobjsense(form) == MinSense ? 1.0 : -1.0

    for res_idx in 1:nb_dual_sols
        if MOI.get(inner, MOI.DualStatus(res_idx)) != MOI.FEASIBLE_POINT
            continue
        end
        solcost = getobjconst(form)
        solconstrs = ConstrId[]
        solvals = Float64[]
        for (id, constr) in getconstrs(form)
            moi_index = getmoiindex(getmoirecord(constr))
            MOI.is_valid(inner, moi_index) || continue
            val = MOI.get(inner, MOI.ConstraintDual(res_idx), moi_index)
            solcost += val * getcurrhs(form, id)
            val = round(val, digits = Coluna.TOL_DIGITS)
            if abs(val) > Coluna.TOL
                push!(solconstrs, id)
                push!(solvals, sense * val)
            end
        end
        varids = VarId[]
        varvals = Float64[]
        activebounds = ActiveBound[]
        for (varid, var) in getvars(form)
            moi_var_index = getmoiindex(getmoirecord(var))
            moi_lower_bound_index = getlowerbound(getmoirecord(var))
            if MOI.is_valid(inner, moi_var_index) && MOI.is_valid(inner, moi_lower_bound_index)
                val = MOI.get(inner, MOI.ConstraintDual(res_idx), moi_lower_bound_index)
                if abs(val) > Coluna.TOL
                    solcost += val * getcurlb(form, varid)
                    push!(varids, varid)
                    push!(varvals, sense * val)
                    push!(activebounds, LOWER)
                end
            end
            moi_upper_bound_index = getupperbound(getmoirecord(var))
            if MOI.is_valid(inner, moi_var_index) && MOI.is_valid(inner, moi_upper_bound_index)
                val = MOI.get(inner, MOI.ConstraintDual(res_idx), moi_upper_bound_index)
                if abs(val) > Coluna.TOL
                    solcost += val * getcurub(form, varid)
                    push!(varids, varid)
                    push!(varvals, sense * val)
                    push!(activebounds, UPPER)
                end
            end
        end
        fixed_obj = 0.0
        for (var_id, fixed_val) in getpartialsol(form)
            cost = getcurcost(form, var_id)
            if abs(cost) > Coluna.TOL
                fixed_obj += cost * fixed_val
            end
        end
        solcost += fixed_obj
        push!(solutions, DualSolution(
            form, solconstrs, solvals, varids, varvals, activebounds, sense*solcost,
            FEASIBLE_SOL
        ))
    end
    return solutions
end

end)  # Core.eval(Coluna.MathProg)

# Use these from model.jl.
lpDualsReset!() = fill!(Coluna.MathProg.lp_duals_stats, 0)
lpDualsPrint() = println("LP duals: fast path ", Coluna.MathProg.lp_duals_stats[1],
                         " calls, original path ", Coluna.MathProg.lp_duals_stats[2], " calls")
