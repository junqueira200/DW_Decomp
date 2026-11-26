# Rank-1 cuts: subset-row cuts with |S| = 3 and multiplier 1/2 (Jepsen et al. 2008),
# added to the Coluna master as non-robust cuts (custom data, see Coluna's tutorial
# docs/src/start/custom_data.jl).
#
#   sum_routes floor(visits of the route to S / 2) * lambda_route <= 1
#
# - Columns carry their route (RouteData); Coluna computes each cut coefficient with
#   computecoeff, same rule as the C++ labeling.
# - Limited memory (Pecin et al. 2017, vertex memory M ⊇ S): walking the route, arriving at
#   a customer outside M resets the counter; a visit to S adds 1, at 2 the coefficient
#   grows by 1 and the counter goes back to 0. The coefficient is <= the full-memory one,
#   so the cut stays valid. M is built from the LP routes so that the violation is the
#   same as with full memory (r1Memory); a cut that the LP escapes later (violated with
#   full memory, not with M) is replaced by one with a larger M.
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
const R1RootOnly        = false  # separate only at the root (rootDone[] = false)
const R1OnlyIfNoCapCuts = false  # separate only in rounds where CVRPSEP found no cut
const R1MaxCutsPerRound = 25     # cuts added per separation round, most violated first
const R1MinViolation    = 0.02   # lhs > 1 + R1MinViolation
const R1DualTol         = 1e-9
const R1LimitedMemory   = true   # false = full memory (M = all customers)
const R1CutMinAge       = 2      # CutCleanup: remove a rank-1 cut after this many calls with dual 0

# Stop rule for the cut rounds (capacity + rank-1), checked at the start of every cut
# callback (capacityCutCallback, model.jl). When it fires, no cut is added, so Coluna ends
# the cut rounds of the node.
#  - time: the last exact pricing call took more than CutStopExactTime seconds (any node);
#  - tailing off (root only): over the last CutTailRounds rounds the master LP value closed
#    less than CutTailGap of the gap (PB - LP) each round on average (VRPSolver: 2%, 3).
const CutStopExactTime = 5.0
const CutTailRounds    = 3
const CutTailGap       = 0.02
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
    M::BitSet                    # memory, S ⊆ M
end

@inline inS(u, S::NTuple{3,Int}) = (u == S[1]) | (u == S[2]) | (u == S[3])

# Full memory: floor(visits to S / 2) (separation: violation of the plain SR cut)
@inline function r1CoefFull(route::Vector{Int}, S::NTuple{3,Int})
    cnt = 0
    for u in route
        cnt += inS(u, S)
    end
    return cnt ÷ 2
end

# Limited memory: the rule of the C++ labeling (state &= keepMask[v], then the S update)
@inline function r1Coef(route::Vector{Int}, S::NTuple{3,Int}, M::BitSet)
    cnt = 0; coef = 0
    for u in route
        u == 0 && continue                       # depot: start (counter 0) or end
        if !(u in M)
            cnt = 0                              # left the memory: forget
        elseif inS(u, S)
            cnt += 1
            cnt == 2 && (coef += 1; cnt = 0)
        end
    end
    return coef
end

# Smallest vertex memory that keeps every route at its full-memory coefficient:
# S + the customers between two visits to S that pair up (counter 1 -> 2).
function r1Memory(S::NTuple{3,Int}, routes)
    M = BitSet(S)
    for r in routes
        open = 0                                 # position of the unpaired visit to S
        for (q, u) in enumerate(r)
            inS(u, S) || continue
            if open == 0
                open = q
            else
                for p in (open + 1):(q - 1)
                    r[p] != 0 && push!(M, r[p])
                end
                open = 0
            end
        end
    end
    return M
end

Coluna.MathProg.computecoeff(v::RouteData, c::R1CutData) = Float64(r1Coef(v.route, c.S, c.M))
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

# Cuts with sigma < 0 in the current pricing round: (S, M, sigma)
const r1Active   = Tuple{NTuple{3,Int},BitSet,Float64}[]
const r1Custs    = Cint[]
const r1Duals    = Float64[]
const r1Ids      = Cint[]
const r1MemStart = Cint[]
const r1MemCusts = Cint[]
const r1Complete = Ref(true)   # false: more cuts with sigma < 0 than the labeling handles

# Stats: pricing calls with cuts, max cuts in one call, separation rounds, cuts added,
# positive duals seen (sign check), calls without exact bound (too many cuts),
# memories enlarged, sum of |M| of the added cuts
const r1Stats = zeros(Int, 8)

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
            σ < -R1DualTol && push!(r1Active, (constr.custom_data.S, constr.custom_data.M, σ))
        end
    else
        for (constrid, σ) in sol
            constrid.custom_family_id == -1 && continue
            constr = MP.getconstr(master, constrid)
            constr.custom_data isa R1CutData || continue
            σ > R1DualTol && (r1Stats[5] += 1)
            σ < -R1DualTol && push!(r1Active, (constr.custom_data.S, constr.custom_data.M, σ))
        end
    end

    empty!(r1Custs); empty!(r1Duals); empty!(r1Ids); empty!(r1MemStart); empty!(r1MemCusts)
    push!(r1MemStart, 0)
    for (k, (S, M, σ)) in enumerate(r1Active)
        push!(r1Custs, S[1], S[2], S[3]); push!(r1Duals, σ); push!(r1Ids, k)
        append!(r1MemCusts, M)
        push!(r1MemStart, length(r1MemCusts))
    end
    loaded = ccall((:setRank1Cuts, LIB), Cint,
                   (Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                   length(r1Duals), r1Custs, r1Duals, r1Ids, r1MemStart, r1MemCusts)
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
    for (S, M, σ) in r1Active
        c = r1Coef(route, S, M)
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

    activeConstrs = Dict{NTuple{3,Int},Any}()      # S -> active constraint
    for (_, constr) in MP.getconstrs(master)
        constr.custom_data isa R1CutData || continue
        MP.iscuractive(master, constr) && (activeConstrs[constr.custom_data.S] = constr)
    end

    # Columns with value > 0: routes and visits per customer
    vals   = Float64[]
    routes = Vector{Vector{Int}}()
    lists  = Vector{Vector{Tuple{Int,Int}}}()      # (customer, visits)
    for (varid, val) in cbdata.orig_sol
        val > 1e-6 || continue
        var = MP.getvar(master, varid)
        var.custom_data isa RouteData || continue
        cnt = Dict{Int,Int}()
        for u in var.custom_data.route
            u != 0 && (cnt[u] = get(cnt, u, 0) + 1)
        end
        push!(vals, val); push!(routes, var.custom_data.route)
        push!(lists, [(u, a) for (u, a) in cnt])
    end
    isempty(vals) && return 0

    lim       = 1.0 + R1MinViolation
    fullM     = BitSet(1:(n - 1))
    memoryOf(S) = R1LimitedMemory ? r1Memory(S, routes) : fullM
    added     = 0
    violEnl   = Float64[]          # true violations (lhs - 1, LP values) of the enlarged cuts
    violNew   = Float64[]          # true violations of the new cuts
    nCands    = 0

    # Active cuts the LP escapes through the memory (violated with full memory, not with M):
    # replaced by a cut with a larger memory; the old one is deactivated.
    if R1LimitedMemory
        for (S, constr) in collect(activeConstrs)
            added >= R1MaxCutsPerRound && break
            M = constr.custom_data.M
            lhsFull = 0.0; lhsM = 0.0
            for c in eachindex(vals)
                lhsFull += vals[c] * r1CoefFull(routes[c], S)
                lhsM    += vals[c] * r1Coef(routes[c], S, M)
            end
            (lhsFull > lim && lhsM <= 1.0 + 1e-6) || continue
            newM = union(M, r1Memory(S, routes))
            lhsNew = 0.0
            for c in eachindex(vals)
                lhsNew += vals[c] * r1Coef(routes[c], S, newM)
            end
            push!(violEnl, lhsNew - 1.0)
            MP.deactivate!(master, constr)
            MOI.submit(model, MOI.UserCut(cbdata),
                       JuMP.ScalarConstraint(JuMP.AffExpr(0.0), MOI.LessThan(1.0)),
                       R1CutData(S, newM))
            added += 1
            r1Stats[7] += 1
            r1Stats[8] += length(newM)
        end
    end

    room = min(R1MaxCutsPerRound - added, R1MaxActive[] - length(activeConstrs))
    if room <= 0
        r1Stats[3] += 1
        r1PrintRound(violNew, violEnl, nCands, length(activeConstrs),
                     added >= R1MaxCutsPerRound ? "per-round limit used by the enlargements"
                                                : "active-cut limit reached")
        return added
    end

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
            (lhs > lim && !haskey(activeConstrs, S)) && push!(cands, (lhs, S))
        end
    end

    sort!(cands, by = first, rev = true)
    nNew = 0
    for (lhs, S) in cands
        nNew >= room && break
        M = memoryOf(S)
        MOI.submit(model, MOI.UserCut(cbdata),
                   JuMP.ScalarConstraint(JuMP.AffExpr(0.0), MOI.LessThan(1.0)),
                   R1CutData(S, M))
        nNew += 1
        r1Stats[8] += length(M)
        push!(violNew, lhs - 1.0)
    end
    nCands = length(cands)
    added += nNew
    r1Stats[3] += 1
    r1Stats[4] += nNew
    r1PrintRound(violNew, violEnl, nCands, length(activeConstrs), "")
    return added
end

# One line per separation round with the true violations (lhs - 1 with the LP values of
# the columns). Coluna's own "avg. viol." line is wrong for these cuts: it adds the
# coefficient of every master column without its value.
function r1PrintRound(violNew, violEnl, nCands, nActive, noRoom::String)
    st(v) = isempty(v) ? "" :
        string(" (avg. viol. ", round(sum(v) / length(v), digits = 3),
               ", max. viol. ", round(maximum(v), digits = 3), ")")
    # noRoom != "": no room for new cuts, the triples were not enumerated
    println("Rank-1 separation: ", length(violNew), " new cuts", st(violNew),
            isempty(noRoom) ? string(" of ", nCands, " violated") : "",
            "; ", length(violEnl), " memories enlarged", st(violEnl),
            "; ", nActive, " active before the round",
            isempty(noRoom) ? "" : string(" (no new cuts: ", noRoom, ")"))
    flush(stdout)                  # progress visible in the log while the run goes on
    return
end

# ---- setup / report --------------------------------------------------------------------

const rootLpHist  = Float64[]   # master LP value at each root cut round
const cutStopStat = zeros(Int, 2) # times fired: [time, tailing off]

function cutStopRule(cbdata)
    lp = Coluna.ColunaBase.getvalue(cbdata.orig_sol)
    if tLastExact[] > CutStopExactTime
        cutStopStat[1] += 1
        println("Cut stop: last exact pricing ", round(tLastExact[], digits = 2), " s > ",
                CutStopExactTime, " s (LP ", round(lp, digits = 4), rootDone[] ? ")" : ", root)")
        flush(stdout)
        return true
    end
    rootDone[] && return false

    push!(rootLpHist, lp)
    k = length(rootLpHist)
    if k > CutTailRounds && isfinite(ubArcFixing[])
        old  = rootLpHist[k - CutTailRounds]
        gap0 = ubArcFixing[] - old
        # average fraction of the gap closed per round over the last CutTailRounds rounds
        closed = gap0 > 1e-9 ? (lp - old) / gap0 / CutTailRounds : 1.0
        if closed < CutTailGap
            cutStopStat[2] += 1
            println("Cut stop: root tailing off, ", round(100 * closed, digits = 2),
                    "% of the gap per round over the last ", CutTailRounds,
                    " rounds (LP ", round(lp, digits = 4), ")")
            flush(stdout)
            return true
        end
    end
    return false
end

function rank1Reset!()
    empty!(rootLpHist)
    fill!(cutStopStat, 0)
    r1SepDualSol[] = nothing
    empty!(r1Active)
    fill!(r1Stats, 0)
    r1Complete[] = true
    R1MaxActive[] = ccall((:getMaxActiveCuts, LIB), Cint, ())
    return
end

function rank1Print()
    println("Cut stop rule: ", cutStopStat[1], " stops by pricing time (> ", CutStopExactTime,
            " s), ", cutStopStat[2], " by root tailing off; root rounds ", length(rootLpHist))
    UseRank1Cuts || return
    println("Rank-1 cuts: ", r1Stats[4], " added in ", r1Stats[3], " separation rounds",
            " (max ", R1MaxCutsPerRound, " per round, ", R1MaxActive[], " active)")
    println("  pricing calls with cuts: ", r1Stats[1], ", max cuts with dual < 0: ", r1Stats[2])
    R1LimitedMemory && println("  limited memory: ", r1Stats[7], " memories enlarged, mean |M| ",
                               round(r1Stats[8] / max(1, r1Stats[4] + r1Stats[7]), digits = 1))
    r1Stats[5] > 0 && println("  WARNING: ", r1Stats[5], " cut duals > 0 (sign convention?)")
    r1Stats[6] > 0 && println("  WARNING: ", r1Stats[6], " pricing calls with too many cuts (no exact bound)")
    return
end
