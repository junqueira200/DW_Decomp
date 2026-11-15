# Cut pool tests (2026-11-15) — not used by the current code

Cut pool for the capacity cuts (one master constraint per customer set S, the callback
skips a known S, CutCleanup reactivates violated inactive cuts). Compared with the
current code (Coluna/CutCleanup.jl, cleanup only, CutMinAge = 5).

| Version | Files | E-n51 (solved) | E-n101 10-min DB |
|---|---|---|---|
| Cleanup only (current) | ../CutCleanup.jl, ../model.jl | 39.3 s, 16 nodes | **805.02** |
| v1: reactivation forces colgen (return true), sticky cuts | CutCleanup_pool_v1.jl + model_pool_v1.jl | 95.1 s, 32 | 803.85 |
| B: return false, sticky cuts | CutCleanup_pool_vB.jl + model_pool_v1.jl | 36.1 s, 14 | 804.26 |
| Non-sticky: return false, CutKeepSticky = false | CutCleanup_pool_nonsticky.jl + model_pool_nonsticky.jl | 34.2 s, 14 | 804.14 |

All pool versions: master LP ~25 ms/solve vs 21.7 ms (cleanup only).

To use one: copy its files over ../CutCleanup.jl and ../model.jl (keep a backup).
The model_*.jl copies are model.jl as it was on 2026-11-15; later changes to model.jl
are not in them.
