cd("/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/SPPRC/Coluna")
using JuMP, BlockDecomposition, Coluna
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using CVRPLIB

include("/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/SPPRC/Coluna/instanceCVRP.jl")
include("/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/SPPRC/Coluna/model.jl")


#global sumWeights = 0.0
CVRPLIB.cache_loaded[] = true; CVRPLIB.instance_id_cache["E-n51-k5"] = 0   # offline: use cached .vrp
global inst::CvrpData = createInstance("E-n51-k5") # E-n76-k7 E-n13-k4 E-n33-k4 

#print("sumWeights: ", sumWeights, "\n\n")
#exit(-1)

global vetReducedCosts::Vector{Float64}       = zeros(inst.n*inst.n)
global vetRoutesOutput::Vector{Cint}          = zeros(Cint, getSizeVetRoutesOutput())
global vetReducedCostsRoutes::Vector{Float64} = zeros(getMaxNumberOfRoutes())
global vetArcMin::Vector{Float64}             = zeros(inst.n*inst.n)    # out of fixArcsSpprc
global removedArcs::BitVector                 = falses(inst.n*inst.n)   # arcs fixed to 0 (whole tree)
# primalBound: known optimum + 1 (E-n76-k7: 682). Change it with the instance, or use Inf.
createModel(inst; primalBound = 522.0)

