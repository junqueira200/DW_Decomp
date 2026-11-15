# Faster column generation helpers for Coluna 0.8.2.
#
# Profile of an E-n101 tree run (2026-11-15): 27% of the time was spent building two
# helper structures at the start of EVERY column generation (each cut round and each
# strong-branching child), Coluna/src/Algorithm/colgen/utils.jl and default.jl:
#   - ReducedCostsCalculationHelper: costs and coefficients of the arc variables
#     (subproblem representatives) in the master constraints. Built twice per column
#     generation: in ColGenContext(...) and again in setup_context!.
#   - SubgradientCalculationHelper: rhs and the same kind of submatrix.
# Their submatrices were built by _submatrix_nz_elems, which reads EVERY master row in
# full (including the entries of all generated columns) and keeps only the arc
# variables: its cost grows with the number of columns.
#
# This file redefines, inside Coluna.Algorithm (same names, same results):
#   (b) _submatrix_nz_elems: reads the columns of the kept variables instead of the rows,
#       so the generated columns are never visited;
#   (a) ReducedCostsCalculationHelper(master) / SubgradientCalculationHelper(master):
#       cached by a key of everything they depend on (the kept variables: active flag and
#       cost; every non-convexity constraint: active, explicit, sense, rhs). Generated
#       columns are not part of the key (the helpers do not depend on them).
# The helpers are only read by Coluna (checked), so a cached one can be shared.
#
# It uses Coluna internals: check it again if Coluna is updated. Julia warns that the
# methods are overwritten.

Core.eval(Coluna.Algorithm, quote

# (b) Same entries as Coluna's version, collected column by column (kept variables only).
function _submatrix_nz_elems(
    form::Formulation,
    keep_constr::Function,
    keep_var::Function,
    m::Function = (form, is_min, constr_id, var_id) -> 1.0
)
    is_min = getobjsense(form) == MinSense
    matrix = getcoefmatrix(form)
    keep_c = Dict{ConstrId,Bool}()
    constr_ids = ConstrId[]
    var_ids = VarId[]
    nz = Float64[]
    for (var_id, var) in getvars(form)
        keep_var(form, var_id, var) || continue
        for (constr_id, coeff) in @view matrix[:, var_id]
            keep = get!(keep_c, constr_id) do
                constr = getconstr(form, constr_id)
                !isnothing(constr) && keep_constr(form, constr_id, constr)
            end
            keep || continue
            push!(constr_ids, constr_id)
            push!(var_ids, var_id)
            push!(nz, m(form, is_min, constr_id, var_id) * coeff)
        end
    end
    return constr_ids, var_ids, nz
end

# (a) Cache of the helpers, keyed by everything they depend on.
const colgen_helper_rch = Dict{UInt,ReducedCostsCalculationHelper}()
const colgen_helper_sh  = Dict{UInt,SubgradientCalculationHelper}()
const colgen_helper_stats = zeros(Int, 4)   # rch hits, rch misses, sh hits, sh misses
const ColGenHelperCacheMax = 500            # entries per cache (emptied when full)

# (c) Row cache (2026-11-15, E-n76-k14 profile: building the helpers was still 22% of
# the run there, many short column generations): the helpers only use the variables
# that represent the arcs (subproblem representatives) and pure master variables.
#  - helper_vars: those variables, collected once (they are created with the model and
#    never change; the master's other variables are the columns, ~10^5 of them);
#  - row cache: for every constraint, its entries on those variables, read once when
#    the constraint is first seen (a constraint's coefficients never change) and kept
#    sorted by variable id.
# The helpers are then built by concatenating the cached rows of the constraints they
# keep, in (constraint, variable) order, with no per-entry dictionary lookups.
is_helper_var_duty(d) = d <= AbstractMasterRepDwSpVar || d <= AbstractOriginMasterVar

const colgen_helper_vars = Tuple{VarId,Variable}[]
const colgen_row_cache   = Dict{ConstrId,Vector{Tuple{VarId,Variable,Float64}}}()

function helper_vars(master)
    if isempty(colgen_helper_vars)
        for (var_id, var) in getvars(master)
            is_helper_var_duty(getduty(var_id)) && push!(colgen_helper_vars, (var_id, var))
        end
        sort!(colgen_helper_vars; by = first)
    end
    return colgen_helper_vars
end

function cached_row(master, constr_id)
    return get!(colgen_row_cache, constr_id) do
        row = Tuple{VarId,Variable,Float64}[]
        for (var_id, coeff) in @view getcoefmatrix(master)[constr_id, :]
            is_helper_var_duty(getduty(var_id)) || continue
            push!(row, (var_id, getvar(master, var_id), coeff))
        end
        sort!(row; by = first)
        row
    end
end

# Same result as _get_costs_and_coeffs(master, keep_var, keep_constr) / _submatrix for
# keep_var predicates that only select helper variables (all the helpers' predicates).
function helper_submatrix(master, keep_constr, keep_var, m = (form, is_min, c, v) -> 1.0)
    is_min = getobjsense(master) == MinSense
    kept = ConstrId[]
    for (constr_id, constr) in getconstrs(master)
        keep_constr(master, constr_id, constr) && push!(kept, constr_id)
    end
    sort!(kept)
    constr_ids = ConstrId[]
    var_ids = VarId[]
    nz = Float64[]
    for constr_id in kept
        for (var_id, var, coeff) in cached_row(master, constr_id)
            keep_var(master, var_id, var) || continue
            push!(constr_ids, constr_id)
            push!(var_ids, var_id)
            push!(nz, m(master, is_min, constr_id, var_id) * coeff)
        end
    end
    return dynamicsparse(
        constr_ids, var_ids, nz, ConstrId(Coluna.MAX_NB_ELEMS), VarId(Coluna.MAX_NB_ELEMS)
    )
end

function helper_costs_and_coeffs(master, keep_var, keep_constr)
    var_ids = VarId[]
    peren_costs = Float64[]
    for (var_id, var) in helper_vars(master)
        if keep_var(master, var_id, var)
            push!(var_ids, var_id)
            push!(peren_costs, getcurcost(master, var_id))
        end
    end
    costs = sparsevec(var_ids, peren_costs, Coluna.MAX_NB_ELEMS)
    return costs, helper_submatrix(master, keep_constr, keep_var)
end

function colgen_helper_key(master)
    h = UInt(0x5bd1e995)
    for (var_id, var) in helper_vars(master)
        h = hash((var_id, iscuractive(master, var), getcurcost(master, var_id)), h)
    end
    for (constr_id, constr) in getconstrs(master)
        getduty(constr_id) <= MasterConvexityConstr && continue
        h = hash((constr_id, iscuractive(master, constr), isexplicit(master, constr),
                  getcursense(master, constr_id), getcurrhs(master, constr_id)), h)
    end
    return h
end

function colgen_helpers_reset!()
    empty!(colgen_helper_rch)
    empty!(colgen_helper_sh)
    empty!(colgen_helper_vars)
    empty!(colgen_row_cache)
    fill!(colgen_helper_stats, 0)
    return
end

# Body of Coluna's ReducedCostsCalculationHelper(master) (colgen/utils.jl), with
# helper_costs_and_coeffs (row cache) instead of _get_costs_and_coeffs.
function build_reduced_costs_helper(master)
    dw_subprob_c, dw_subprob_A = helper_costs_and_coeffs(
        master,
        (form, var_id, var) -> getduty(var_id) <= AbstractMasterRepDwSpVar && iscuractive(form, var),
        (form, constr_id, constr) -> !(getduty(constr_id) <= MasterConvexityConstr) && iscuractive(form, constr)
    )
    master_c, master_A = helper_costs_and_coeffs(
        master,
        (form, var_id, var) -> getduty(var_id) <= AbstractOriginMasterVar && iscuractive(form, var),
        (form, constr_id, constr) -> !(getduty(constr_id) <= MasterConvexityConstr)
    )
    return ReducedCostsCalculationHelper(dw_subprob_c, dw_subprob_A, master_c, master_A)
end

# Body of Coluna's SubgradientCalculationHelper(master) (colgen/utils.jl), with
# helper_submatrix (row cache) instead of _submatrix.
function build_subgradient_helper(master)
    m_rhs = (master, is_min, constr_id) -> begin
        constr_sense = getcursense(master, constr_id)
        if is_min
            return constr_sense == Less ? -1.0 : 1.0
        else
            return constr_sense == Greater ? -1.0 : 1.0
        end
    end
    m_submatrix = (master, is_min, constr_id, var_id) -> begin
        m_rhs(master, is_min, constr_id)
    end

    constr_ids = ConstrId[]
    constr_rhs = Float64[]
    constr_rhs_dual = Float64[]

    is_min = getobjsense(master) == MinSense
    for (constr_id, constr) in getconstrs(master)
        if !(getduty(constr_id) <= MasterConvexityConstr) &&
           iscuractive(master, constr) && isexplicit(master, constr)
            push!(constr_ids, constr_id)
            push!(constr_rhs, m_rhs(master, is_min, constr_id) * getcurrhs(master, constr_id))
            push!(constr_rhs_dual, getcurrhs(master, constr_id))
        end
    end

    a = sparsevec(constr_ids, constr_rhs, Coluna.MAX_NB_ELEMS)
    a_dual = sparsevec(constr_ids, constr_rhs_dual, Coluna.MAX_NB_ELEMS)
    A = helper_submatrix(
        master,
        (form, constr_id, constr) -> !(getduty(constr_id) <= MasterConvexityConstr) && iscuractive(form, constr),
        (form, var_id, var) -> getduty(var_id) <= MasterPureVar || getduty(var_id) <= MasterRepPricingVar,
        m_submatrix
    )
    return SubgradientCalculationHelper(a, a_dual, A)
end

function ReducedCostsCalculationHelper(master)
    key = colgen_helper_key(master)
    h = get(colgen_helper_rch, key, nothing)
    if !isnothing(h)
        colgen_helper_stats[1] += 1
        return h
    end
    colgen_helper_stats[2] += 1
    length(colgen_helper_rch) >= ColGenHelperCacheMax && empty!(colgen_helper_rch)
    h = build_reduced_costs_helper(master)
    colgen_helper_rch[key] = h
    return h
end

function SubgradientCalculationHelper(master)
    key = colgen_helper_key(master)
    h = get(colgen_helper_sh, key, nothing)
    if !isnothing(h)
        colgen_helper_stats[3] += 1
        return h
    end
    colgen_helper_stats[4] += 1
    length(colgen_helper_sh) >= ColGenHelperCacheMax && empty!(colgen_helper_sh)
    h = build_subgradient_helper(master)
    colgen_helper_sh[key] = h
    return h
end

end)  # Core.eval(Coluna.Algorithm)

# Use these from model.jl.
colgenHelpersReset!() = Coluna.Algorithm.colgen_helpers_reset!()
function colgenHelpersPrint()
    s = Coluna.Algorithm.colgen_helper_stats
    println("ColGen helpers cache: reduced costs ", s[1], " hits / ", s[2], " builds; ",
            "subgradient ", s[3], " hits / ", s[4], " builds")
end
