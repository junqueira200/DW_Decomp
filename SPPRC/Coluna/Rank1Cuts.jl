# Rank-1 cuts: subset-row cuts with |S| = 3 and multiplier 1/2 (Jepsen et al. 2008),
# added to the Coluna master as non-robust cuts (custom data, see Coluna's tutorial
# docs/src/start/custom_data.jl).
#
#   sum_routes floor(visits of the route to S / 2) * lambda_route <= 1
#
# - Columns carry their route (RouteData); Coluna computes each cut coefficient with
#   computecoeff (same rule as the C++ labeling: count visits to S, divide by 2, round down;
#   full memory, no limited memory yet).
# - Pricing: the cut duals sigma <= 0 are passed to the C++ library (setRank1Cuts) before
#   every solveSpprc call; the labeling charges -sigma each time a route completes 2
#   visits to S. The C++ library does not return the route reduced costs, so the
#   penalties are added here, in submitRoutesToColuna.
# - Separation: in the cut callback, after CVRPSEP. All triples of customers are checked
#   against the master solution; the most violated ones are added.
#
# Duals: with smoothing stabilization the pricing must use the smoothed duals (the arc
# reduced costs given to the callback are computed with them). Coluna stores the real LP
# duals in the master constraints (getcurincval), so the smoothed dual solution is taken
# from get_stab_dual_sol (wrapped below) to keep arcs and cuts consistent. Coluna computes
# the real reduced cost of the columns itself (_nonrobust_cuts_contrib).

const UseRank1Cuts      = true
const R1RootOnly        = true   # separate only at the root (rootDone[] = false)
const R1OnlyIfNoCapCuts = false  # separate only in rounds where CVRPSEP found no cut
const R1MaxCutsPerRound = 10     # cuts added per separation round, most violated first
const R1MinViolation    = 0.02   # lhs > 1 + R1MinViolation
const R1DualTol         = 1e-9
# At most this many rank-1 cuts are active in the master: the C++ labeling keeps one bit
# per cut with sigma < 0 (MaxActiveCuts in LabelingConstants.h, read in rank1Reset!).
const R1MaxActive       = Ref(0)

struct RouteData <: BD.AbstractCustomVarData
    route::Vector{Int}           # 0 c1 ... ck 0
end
# Coluna's pool compares custom data with == (same column = same route)
Base.:(==)(a::RouteData, b::RouteData) = a.route == b.route
Base.hash(a::RouteData, h::UInt) = hash(a.route, hash(:RouteData, h))

struct R1CutData <: BD.AbstractCustomConstrData
    S::NTuple{3,Int}             # customers, sorted
end

@inline function r1Coef(route::Vector{Int}, S::NTuple{3,Int})
    cnt = 0
    for u in route
        cnt += (u == S[1]) | (u == S[2]) | (u == S[3])
    end
    return cnt ÷ 2
end

Coluna.MathProg.computecoeff(v::RouteData, c::R1CutData) = Float64(r1Coef(v.route, c.S))
# Coluna's MOI.submit of a user cut calls computecoeff for every column with custom data,
# also for robust cuts (capacity cuts, custom data = nothing): no non-robust part.
Coluna.MathProg.computecoeff(::RouteData, ::Nothing) = 0.0

# ---- duals used by the pricing ---------------------------------------------------------

const r1SepDualSol = Ref{Any}(nothing)     # dual solution given to the last pricing round

function Coluna.ColGen.get_stab_dual_sol(stab::Coluna.Algorithm.ColGenStab,
                                         phase::Coluna.ColGen.AbstractColGenPhase,
                                         mast_dual_sol::Coluna.MathProg.DualSolution)
    sol = invoke(Coluna.ColGen.get_stab_dual_sol,
                 Tuple{Coluna.Algorithm.ColGenStab, Any, Any}, stab, phase, mast_dual_sol)
    r1SepDualSol[] = sol
    return sol
end

function Coluna.ColGen.get_stab_dual_sol(stab::Coluna.Algorithm.NoColGenStab,
                                         phase::Coluna.ColGen.AbstractColGenPhase,
                                         mast_dual_sol::Coluna.MathProg.DualSolution)
    r1SepDualSol[] = mast_dual_sol
    return mast_dual_sol
end

# Cuts with sigma < 0 in the current pricing round: (S, sigma)
const r1Active   = Tuple{NTuple{3,Int},Float64}[]
const r1Custs    = Cint[]
const r1Duals    = Float64[]
const r1Ids      = Cint[]
const r1Complete = Ref(true)   # false: more cuts with sigma < 0 than the labeling handles

# Stats: pricing calls with cuts, max cuts in one call, separation rounds, cuts added,
# positive duals seen (sign check), calls without exact bound (too many cuts)
const r1Stats = zeros(Int, 6)

# Reads the cut duals and loads them in the C++ library. Call before every solveSpprc.
function r1LoadPricingCuts!(cbdata)
    empty!(r1Active)
    UseRank1Cuts || return
    MP     = Coluna.MathProg
    master = cbdata.form.parent_formulation
    sol    = r1SepDualSol[]

    if sol === nothing
        for (_, constr) in MP.getconstrs(master)
            constr.custom_data isa R1CutData || continue
            σ = MP.getcurincval(master, constr)
            σ > R1DualTol && (r1Stats[5] += 1)
            σ < -R1DualTol && push!(r1Active, (constr.custom_data.S, σ))
        end
    else
        for (constrid, σ) in sol
            constrid.custom_family_id == -1 && continue
            constr = MP.getconstr(master, constrid)
            constr.custom_data isa R1CutData || continue
            σ > R1DualTol && (r1Stats[5] += 1)
            σ < -R1DualTol && push!(r1Active, (constr.custom_data.S, σ))
        end
    end

    empty!(r1Custs); empty!(r1Duals); empty!(r1Ids)
    for (k, (S, σ)) in enumerate(r1Active)
        push!(r1Custs, S[1], S[2], S[3]); push!(r1Duals, σ); push!(r1Ids, k)
    end
    loaded = ccall((:setRank1Cuts, LIB), Cint, (Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}),
                   length(r1Duals), r1Custs, r1Duals, r1Ids)
    r1Complete[] = loaded >= 0
    r1Complete[] || (r1Stats[6] += 1)
    if !isempty(r1Active)
        r1Stats[1] += 1
        r1Stats[2] = max(r1Stats[2], length(r1Active))
    end
    return
end

# Reduced-cost contribution of the cuts to a route: -sum sigma * coef >= 0
function r1RoutePenalty(route::Vector{Int})
    pen = 0.0
    for (S, σ) in r1Active
        c = r1Coef(route, S)
        c > 0 && (pen -= σ * c)
    end
    return pen
end

# ---- separation ------------------------------------------------------------------------

# Adds up to R1MaxCutsPerRound violated cuts; returns the number added.
function r1Separate(cbdata, model, n::Int)
    UseRank1Cuts || return 0
    (R1RootOnly && rootDone[]) && return 0
    MP     = Coluna.MathProg
    master = cbdata.form

    active = Set{NTuple{3,Int}}()
    for (_, constr) in MP.getconstrs(master)
        constr.custom_data isa R1CutData || continue
        MP.iscuractive(master, constr) && push!(active, constr.custom_data.S)
    end
    room = min(R1MaxCutsPerRound, R1MaxActive[] - length(active))
    room <= 0 && return 0

    # Columns with value > 0: visits per customer
    vals  = Float64[]
    lists = Vector{Vector{Tuple{Int,Int}}}()       # (customer, visits)
    for (varid, val) in cbdata.orig_sol
        val > 1e-6 || continue
        var = MP.getvar(master, varid)
        var.custom_data isa RouteData || continue
        cnt = Dict{Int,Int}()
        for u in var.custom_data.route
            u != 0 && (cnt[u] = get(cnt, u, 0) + 1)
        end
        push!(vals, val); push!(lists, [(u, a) for (u, a) in cnt])
    end
    isempty(vals) && return 0

    # Upper bound of the lhs of S = {i,j,k}: W(i,j) + W(i,k) + W(j,k) + D(i) + D(j) + D(k)
    # with W(i,j) = sum val*a_i*a_j and D(i) = sum val*floor(a_i/2)
    # (floor((a+b+c)/2) <= ab + ac + bc + floor(a/2) + floor(b/2) + floor(c/2)).
    W = zeros(n, n); D = zeros(n)
    V = zeros(Int, length(vals), n)
    for (c, lst) in enumerate(lists)
        v = vals[c]
        for (i, a) in lst
            V[c, i] = a
            D[i] += v * (a ÷ 2)
            for (j, b) in lst
                i < j && (W[i, j] += v * a * b; W[j, i] = W[i, j])
            end
        end
    end

    lim   = 1.0 + R1MinViolation
    cands = Tuple{Float64,NTuple{3,Int}}[]
    for i in 1:(n - 1), j in (i + 1):(n - 1)
        wij = W[i, j] + D[i] + D[j]
        for k in (j + 1):(n - 1)
            wij + W[i, k] + W[j, k] + D[k] > lim || continue
            lhs = 0.0
            for c in eachindex(vals)
                s = V[c, i] + V[c, j] + V[c, k]
                s >= 2 && (lhs += vals[c] * (s ÷ 2))
            end
            S = (i, j, k)
            (lhs > lim && !(S in active)) && push!(cands, (lhs, S))
        end
    end

    sort!(cands, by = first, rev = true)
    added = 0
    for (lhs, S) in cands
        added >= room && break
        MOI.submit(model, MOI.UserCut(cbdata),
                   JuMP.ScalarConstraint(JuMP.AffExpr(0.0), MOI.LessThan(1.0)),
                   R1CutData(S))
        added += 1
    end
    r1Stats[3] += 1
    r1Stats[4] += added
    printInfo && println("rank-1 cuts: ", length(cands), " violated, ", added, " added")
    return added
end

# ---- setup / report --------------------------------------------------------------------

function rank1Reset!()
    r1SepDualSol[] = nothing
    empty!(r1Active)
    fill!(r1Stats, 0)
    r1Complete[] = true
    R1MaxActive[] = ccall((:getMaxActiveCuts, LIB), Cint, ())
    return
end

function rank1Print()
    UseRank1Cuts || return
    println("Rank-1 cuts: ", r1Stats[4], " added in ", r1Stats[3], " separation rounds",
            " (max ", R1MaxCutsPerRound, " per round, ", R1MaxActive[], " active)")
    println("  pricing calls with cuts: ", r1Stats[1], ", max cuts with dual < 0: ", r1Stats[2])
    r1Stats[5] > 0 && println("  WARNING: ", r1Stats[5], " cut duals > 0 (sign convention?)")
    r1Stats[6] > 0 && println("  WARNING: ", r1Stats[6], " pricing calls with too many cuts (no exact bound)")
    return
end
