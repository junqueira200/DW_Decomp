# PGO training run (called by pgo/pgo_train.sh) with Coluna/model.jl.
#
# usage: julia [-J sysimage] train.jl <instance> <primal bound> [tree seconds]
#   without [tree seconds]: root node only (no strong branching, one node);
#   with [tree seconds]:    the normal tree search (strong branching), stopped after
#                           that many seconds (TimeLimitSec).
#
# Run it from a folder whose ../build/libSPPRC.dylib is the library to profile
# (Coluna/instanceCVRP.jl: const LIB = "../build/libSPPRC.dylib").

using JuMP, BlockDecomposition, Coluna
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using CVRPLIB

const colunaDir   = normpath(joinpath(@__DIR__, "..", "Coluna"))
const instName    = ARGS[1]
const primalBound = parse(Float64, ARGS[2])
const treeSeconds = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0

include(joinpath(colunaDir, "instanceCVRP.jl"))

# model.jl changed for the training run. If a pattern is not found, model.jl was
# changed: update the pattern here.
let src = read(joinpath(colunaDir, "model.jl"), String)
    changes = treeSeconds > 0 ?
        ["const TimeLimitSec = 4.0*60.0" => "const TimeLimitSec = $(treeSeconds)"] :
        ["global UseStrongBranching::Bool = true" => "global UseStrongBranching::Bool = false",
         "maxnumnodes = 10000000"                 => "maxnumnodes = 1"]
    for (a, b) in changes
        occursin(a, src) || error("train.jl: \"$a\" not found in model.jl")
        src = replace(src, a => b)
    end
    include_string(Main, src, joinpath(colunaDir, "model.jl"))
end

CVRPLIB.cache_loaded[] = true; CVRPLIB.instance_id_cache[instName] = 0   # offline: cached .vrp
global inst::CvrpData = createInstance(instName)

global vetReducedCosts::Vector{Float64}       = zeros(inst.n*inst.n)
global vetRoutesOutput::Vector{Cint}          = zeros(Cint, getSizeVetRoutesOutput())
global vetReducedCostsRoutes::Vector{Float64} = zeros(getMaxNumberOfRoutes())
global vetArcMin::Vector{Float64}             = zeros(inst.n*inst.n)    # out of fixArcsSpprc
global removedArcs::BitVector                 = falses(inst.n*inst.n)   # arcs fixed to 0
createModel(inst; primalBound = primalBound)
