# Rank-1 cuts and their effect on the pricing problem

Notes written on 2026-11-16. Context: E-n76-k7 has a root bound of 669.02 with the
current code, against an optimum of 682 (gap 12.98, 1.9%). After 4 minutes of tree
search the bound is about 672.9 (1.3% gap). Most of the gap comes from the root bound;
rank-1 cuts are what state-of-the-art CVRP solvers use to close it.

## What rank-1 cuts are

They are Chvátal–Gomory cuts of rank 1, built on the set-partitioning rows of the master
(one row per customer: "each customer is visited exactly once").

**Subset-row cuts (SRC, Jepsen et al., 2008)** are the most used kind. For a subset S of
three customers and multiplier ½:

    Σ_routes ⌊½ · (number of customers of S the route visits)⌋ · λ_route ≤ 1

- A route visiting 2 or 3 customers of S has coefficient 1.
- A route visiting 0 or 1 of them has coefficient 0.
- This cuts fractional LP solutions where several routes share pairs of the three
  customers, which capacity cuts cannot cut.

The general rank-1 form uses larger subsets (up to 5) and other multiplier vectors
(Pecin et al.).

**Strength:** Jepsen et al. report that these cuts raise the lower bound significantly
and often prove optimality at the root. Together with limited memory, they are the core
of the BaPCod / VRPSolver results.

**They are non-robust.** A route's coefficient depends on how many customers of S it
visits, not on the arcs it uses. So, unlike the capacity cuts, their duals cannot be
folded into the arc reduced costs: the pricing problem itself has to change.

## How the pricing (labeling) changes

For each active cut c with dual σ_c ≤ 0 (½-multiplier subset-row case):

1. **New label state:** each label carries, for each cut, a counter s_c ∈ {0, 1}
   (customers of S_c visited so far, mod 2).
2. **Extension to customer j:** for every cut with j ∈ S_c, s_c increases by 1. When it
   reaches 2, it goes back to 0 and the label pays the penalty: its reduced cost
   increases by |σ_c|.
3. **Dominance:** besides reduced cost, load and ng-sets, L₁ dominates L₂ only if

       rc(L₁) + Σ_{c : s₁,c > s₂,c} |σ_c| ≤ rc(L₂)

   L₁ is charged the penalty for each cut where it is one step closer to paying than L₂.
   This is Jepsen et al.'s modified dominance rule, which uses the step-like cost
   structure.
4. **Limited memory (Pecin, Pessoa, Poggi, Uchoa and Santos, 2017).** This is what makes
   the cuts workable in practice:
   - each cut has a memory set M_c ⊇ S_c (customers or arcs);
   - when a label moves outside M_c, its counter s_c resets to 0;
   - the cut becomes slightly weaker but stays valid;
   - many more labels end up with equal counters, so dominance stays strong;
   - the memory sets are chosen from the current LP solution and enlarged when needed.

   Without limited memory, each cut adds a real dimension to dominance, and labeling
   with many cuts blows up.

## What it would mean for this code

### C++ labeling (libSPPRC)

- **Label:** a small array or bitset of counters for the active cuts.
- **Extension:** update the counters and add the penalties.
- **Dominance:** in `dominanceIntraBucketSlow` and in the check against extended labels,
  add the cut-penalty term. This makes dominance more expensive and weaker; that is the
  price of the stronger bound.
- **Existing parts:**
  - **heuristic stages:** they can ignore the cuts in dominance, the same way they
    already ignore ng-sets with `useNg = false`;
  - **completion bound (`matDp`) and arc fixing:** they stay valid if they ignore the
    cuts, since the penalties are never negative. They only get somewhat weaker.
- **C API:** each pricing call would also pass the active cuts (subsets, memory sets,
  multipliers) and their duals.

### Coluna side

- **Non-robust cuts are supported:** they need custom data on the cuts *and* on every
  column (BlockDecomposition custom data plus a `computecoeff` method). Coluna's
  tutorial shows this with bin packing.
- **Columns:** they would need to carry their route, or at least their customer visit
  counts, so `computecoeff` can compute each cut's coefficient.
- **Pricing callback:** it reads the cut duals from the master constraints (as in the
  tutorial) and passes them to C++.
- **Separation:** add a subset-row separator (enumerating triples is cheap for n ≈ 100)
  to the cut callback, and submit the cuts with their custom data.
- **Performance risk:** for non-robust constraints, Coluna computes coefficients by
  looping over all variables and constraints (`_computenonrobustmembers`). That is the
  same kind of loop patched in `ColGenHelpers.jl`, `ActiveSets.jl` and `CutInsert.jl`,
  and it would probably need a patch too.
- **Maturity:** Coluna's README lists cut generation (robust and non-robust) as working
  but still needing more testing and performance review.

### Effort and a first step

A large project, mostly in the C++ labeling: new label state, a new dominance rule,
limited-memory management, and the C API. On the Julia side: custom data, a separator,
and probably one or two Coluna patches.

A reasonable first step: plain subset-row cuts on three customers, with limited memory
on customers only, separated at the root only, to measure how much the root gap closes
before building more.

A reference implementation exists in the bucket-graph library
(<https://github.com/spoorendonk/bucket-graph-spprc>): `include/bgspprc/r1c.h`
(about 1,200 lines) implements limited-memory rank-1 cuts as a labeling resource.

## Sources

- [Subset-Row Inequalities Applied to the VRPTW (Jepsen, Petersen, Spoorendonk, Pisinger, Operations Research 56(2), 2008), IDEAS/RePEc](https://ideas.repec.org/a/inm/oropre/v56y2008i2p497-511.html)
- [Subset-row inequalities applied to the VRPTW (University of Copenhagen record)](https://test01.ku.dk/:obvius/pureproxy/3974268/en/publications/subset-row-inequalities-applied-to-the-vehicle-routing-problem-wi)
- [Limited memory rank-1 cuts for vehicle routing problems (Pecin et al., Operations Research Letters 45(3), 2017), EconBiz](https://econbiz.de/Record/limited-memory-rank-1-cuts-for-vehicle-routing-problems-pecin-diego/10011719279)
- [Limited memory Rank-1 Cuts for Vehicle Routing Problems (Bohrium summary)](https://www.bohrium.com/paper-details/limited-memory-rank-1-cuts-for-vehicle-routing-problems/811058978619392000-11870)
- [RouteOpt documentation: Rank-1 Cuts](https://zhengzhong-you.github.io/RouteOpt-Docs/rank1_cuts/r1c_index.html)
- [Enhancements of Fragment Based Algorithms for Vehicle Routing Problems (arXiv 2411.13151)](https://arxiv.org/pdf/2411.13151)
- [Sadykov, VeRoLog 2017 slides](https://www.math.u-bordeaux.fr/~rsadykov/slides/Sadykov_Verolog17slides.pdf)
- [Pessoa, Sadykov, Uchoa, EJOR 2018 preprint](https://www.math.u-bordeaux.fr/~rsadykov/papers/PessoaSadykovUchoa_EJOR18.pdf)
- [Mainz discussion paper (improved dominance with subset-row cuts)](https://logistik.bwl.uni-mainz.de/files/2025/02/LM-2023-04.R1.pdf)
- [Coluna.jl: Custom Variables and Cuts (non-robust cuts)](https://atoptima.github.io/Coluna.jl/dev/start/custom_data)
- [Coluna.jl on GitHub](https://www.Github.com/atoptima/Coluna.jl)
