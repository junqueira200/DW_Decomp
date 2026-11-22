# Extra sections in Coluna's timer table (Coluna._to, printed at the end of the run).
#
# Coluna 0.8.2 only times "Pricing callback" and "SolveLpForm"; on an E-n101-k8 root
# about 41% of the time was not measured. This file adds:
#  - Coluna's main steps: column generation (each iteration split into master LP,
#    reduced costs, pricing, column insertion, ...), cut separation, cut cleanup / arc
#    fixing, primal heuristics, node records and branching;
#  - the macro @ctime "name" expr, used in model.jl inside the pricing and cut callbacks.
#
# How Coluna's functions are timed: each hook is a NEW method that is slightly more
# specific than Coluna's one (one argument that Coluna leaves untyped gets a type). It
# times the call and runs Coluna's own method through `invoke`. Nothing of Coluna is
# copied or overwritten, so Julia gives no "method overwritten" warning. If a hook does
# not match (e.g. another argument type), Coluna's method runs untimed: never an error.
# The hooks use Coluna internals (function names and argument lists): check them again
# if Coluna is updated.
#
# UseTimers = false removes every timer (the hooks are not defined and @ctime expands to
# the bare expression). The cost of a timer is well under a microsecond.

const UseTimers = true

macro ctime(name, ex)
    UseTimers || return esc(ex)
    return esc(:(Coluna.TO.@timeit Coluna._to $name $ex))
end

# Julia's JIT compilation during optimize! (code that is not in the sysimage, or that
# must be compiled again because a patch file / hook changed a method it calls). It is
# included in the times of the table: compiling happens the first time a function runs,
# inside whatever section called it.
const compileAtStart = Ref(0.0)
function timersStart!()
    Base.cumulative_compile_timing(true)
    compileAtStart[] = Base.cumulative_compile_time_ns()[1] / 1e9
end
function timersPrint()
    t = Base.cumulative_compile_time_ns()[1] / 1e9 - compileAtStart[]
    println("JIT compilation during optimize!: ", round(t, digits = 2),
            " s (already included in the times of the table)")
end

if UseTimers

const CA   = Coluna.Algorithm
const CG   = Coluna.ColGen
const CMP  = Coluna.MathProg
const PCtx = Coluna.Algorithm.ColGenPrinterContext   # column generation context (print = true)
const CCtx = Coluna.Algorithm.ColCutGenContext       # conquer context

# ---------------------------------------------------------------- tree and conquer

# Initialization of the storage units (once, before the root).
function CA.initialize_storage_units!(reform::CMP.Reformulation{CMP.DwMaster}, algo::CA.AbstractOptimizationAlgorithm)
    @ctime "Storage init" invoke(CA.initialize_storage_units!,
                                 Tuple{CMP.Reformulation, CA.AbstractOptimizationAlgorithm},
                                 reform, algo)
end

# Tree search: creation of the search space (once).
function Coluna.TreeSearch.new_space(T::Type{CA.BaBSearchSpace}, algo::CA.TreeSearchAlgorithm,
                                     reform::CMP.Reformulation, input::CA.OptimizationState)
    @ctime "Tree setup" invoke(Coluna.TreeSearch.new_space,
                               Tuple{Type{CA.BaBSearchSpace}, CA.TreeSearchAlgorithm, CMP.Reformulation, Any},
                               T, algo, reform, input)
end

# One node of the tree: restore of the node's state (node_change!), conquer, after
# conquer, branching. Its time outside "Conquer" and "Branching" is the node change.
function Coluna.TreeSearch.children(space::CA.BaBSearchSpace, current::CA.Node, env::Coluna.Env)
    @ctime "Node" invoke(Coluna.TreeSearch.children,
                         Tuple{CA.AbstractColunaSearchSpace, Coluna.TreeSearch.AbstractNode, Any},
                         space, current, env)
end

function CA.after_conquer!(space::CA.BaBSearchSpace, current::CA.Node, conquer_output)
    @ctime "After conquer" invoke(CA.after_conquer!, Tuple{CA.BaBSearchSpace, Any, Any},
                                  space, current, conquer_output)
end

# Conquer of one node (column generation, cut rounds, heuristics).
function CA.run_colcutgen_conquer!(ctx::CCtx, env::Coluna.Env, reform, input)
    @ctime "Conquer" invoke(CA.run_colcutgen_conquer!, Tuple{CCtx, Any, Any, Any},
                            ctx, env, reform, input)
end


# Column generation of one cut round (also builds the colgen context and helpers).
function CA.run_colgen!(ctx::CCtx, env::Coluna.Env, reform, input, conquer_output)
    @ctime "Column generation" invoke(CA.run_colgen!, Tuple{CCtx, Any, Any, Any, Any},
                                      ctx, env, reform, input, conquer_output)
end

# Algorithm run before each cut round: CutCleanup (cut pool cleanup, column cleanup and
# node arc fixing, see CutCleanup.jl).
function CA.run_before_cutgen_user_algo!(ctx::CCtx, algo, env::Coluna.Env, reform, conquer_output)
    @ctime "Cut cleanup + arc fixing" invoke(CA.run_before_cutgen_user_algo!,
                                             Tuple{CCtx, Any, Any, Any, Any},
                                             ctx, algo, env, reform, conquer_output)
end

# One cut round: separation (capacityCutCallback) and insertion of the cuts in the master.
function CA.run_cutgen!(ctx::CCtx, env::Coluna.Env, reform, sol)
    @ctime "Cut separation round" invoke(CA.run_cutgen!, Tuple{CCtx, Any, Any, Any},
                                         ctx, env, reform, sol)
end

# Primal heuristics (restricted master IP at the root).
function CA.run_heuristics!(ctx::CCtx, heuristics, env::Coluna.Env, reform, input, conquer_output)
    @ctime "Primal heuristics" invoke(CA.run_heuristics!, Tuple{CCtx, Any, Any, Any, Any, Any},
                                      ctx, heuristics, env, reform, input, conquer_output)
end

# Records of the formulations (node state, for the tree and strong branching).
function CA.create_records(reform::CMP.Reformulation{CMP.DwMaster})
    @ctime "Node records" invoke(CA.create_records, Tuple{CMP.Reformulation}, reform)
end

# Branching (strong branching includes the evaluation of the candidates).
function CA.run!(algo::CA.StrongBranching, env::Coluna.Env, reform::CMP.Reformulation,
                 input::Coluna.Branching.AbstractDivideInput)
    @ctime "Branching" invoke(CA.run!, Tuple{Coluna.AlgoAPI.AbstractDivideAlgorithm, Coluna.Env,
                                             CMP.Reformulation, Coluna.Branching.AbstractDivideInput},
                              algo, env, reform, input)
end
function CA.run!(algo::CA.ClassicBranching, env::Coluna.Env, reform::CMP.Reformulation,
                 input::Coluna.Branching.AbstractDivideInput)
    @ctime "Branching" invoke(CA.run!, Tuple{Coluna.AlgoAPI.AbstractDivideAlgorithm, Coluna.Env,
                                             CMP.Reformulation, Coluna.Branching.AbstractDivideInput},
                              algo, env, reform, input)
end

# ---------------------------------------------------------------- column generation

function CG.run!(ctx::PCtx, env, ip_primal_sol; iter = 1)
    @ctime "colgen loop" invoke(CG.run!, Tuple{Any, Any, Any}, ctx, env, ip_primal_sol; iter = iter)
end

function CG.run_colgen_iteration!(ctx::PCtx, phase, stage, env, ip_primal_sol, stab)
    @ctime "iteration" invoke(CG.run_colgen_iteration!, Tuple{Any, Any, Any, Any, Any, Any},
                              ctx, phase, stage, env, ip_primal_sol, stab)
end

# Master LP (contains Coluna's "SolveLpForm": the rest is Coluna's handling around it).
function CG.optimize_master_lp_problem!(master::CMP.Formulation, ctx::PCtx, env)
    @ctime "master LP" invoke(CG.optimize_master_lp_problem!, Tuple{Any, PCtx, Any},
                              master, ctx, env)
end

# Integrality check of the master LP solution (and essential cuts, if any).
function CG.check_primal_ip_feasibility!(sol, ctx::PCtx, phase, env::Coluna.Env)
    @ctime "IP feasibility check" invoke(CG.check_primal_ip_feasibility!, Tuple{Any, PCtx, Any, Any},
                                         sol, ctx, phase, env)
end

function CG.update_stabilization_after_master_optim!(stab::CA.ColGenStab, phase::CG.AbstractColGenPhase, dual_sol)
    @ctime "stabilization" invoke(CG.update_stabilization_after_master_optim!,
                                  Tuple{CA.ColGenStab, Any, Any}, stab, phase, dual_sol)
end

function CG.update_master_constrs_dual_vals!(ctx::PCtx, dual_sol::CMP.DualSolution)
    @ctime "store duals" invoke(CG.update_master_constrs_dual_vals!, Tuple{PCtx, Any},
                                ctx, dual_sol)
end

function CG.update_reduced_costs!(ctx::PCtx, phase::CG.AbstractColGenPhase, red_costs)
    @ctime "reduced costs (Coluna)" invoke(CG.update_reduced_costs!, Tuple{PCtx, Any, Any},
                                           ctx, phase, red_costs)
end

function CG.update_sp_vars_red_costs!(ctx::PCtx, sp::CMP.Formulation{CMP.DwSp}, red_costs::AbstractVector)
    @ctime "reduced costs (Coluna)" invoke(CG.update_sp_vars_red_costs!,
                                           Tuple{PCtx, CMP.Formulation{CMP.DwSp}, Any},
                                           ctx, sp, red_costs)
end

# What is left of "iteration" outside its sections is mostly Coluna's computation of the
# reduced costs of all subproblem variables, c - A' * duals (ColGen.jl, no function to
# hook), once per iteration and again when stabilization changes the duals.

# Pricing (contains Coluna's "Pricing callback").
function CG.optimize_pricing_problem!(ctx::PCtx, sp::CMP.Formulation{CMP.DwSp}, env::Coluna.Env,
                                      optimizer, master_dual_sol, stab_changes_mast_dual_sol)
    @ctime "pricing" invoke(CG.optimize_pricing_problem!,
                            Tuple{PCtx, CMP.Formulation{CMP.DwSp}, Any, Any, Any, Any},
                            ctx, sp, env, optimizer, master_dual_sol, stab_changes_mast_dual_sol)
end

function CG.push_in_set!(ctx::PCtx, set::CA.ColumnsSet, col)
    @ctime "push columns" invoke(CG.push_in_set!, Tuple{PCtx, Any, Any}, ctx, set, col)
end

function CG.compute_dual_bound(ctx::PCtx, phase::CG.AbstractColGenPhase, sp_dbs, generated_columns, master_dual_sol)
    @ctime "dual bound" invoke(CG.compute_dual_bound, Tuple{PCtx, Any, Any, Any, Any},
                               ctx, phase, sp_dbs, generated_columns, master_dual_sol)
end

function CG.update_stabilization_after_pricing_optim!(stab, ctx::PCtx, generated_columns,
                                                      master::CMP.Formulation, pseudo_db, smooth_dual_sol)
    @ctime "stabilization" invoke(CG.update_stabilization_after_pricing_optim!,
                                  Tuple{Any, PCtx, Any, Any, Any, Any},
                                  stab, ctx, generated_columns, master, pseudo_db, smooth_dual_sol)
end

function CG.insert_columns!(ctx::PCtx, phase::CG.AbstractColGenPhase, columns)
    @ctime "insert columns" invoke(CG.insert_columns!, Tuple{PCtx, Any, Any}, ctx, phase, columns)
end

end # if UseTimers
