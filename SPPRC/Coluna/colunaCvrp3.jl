using JuMP, BlockDecomposition, Coluna, Gurobi
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using CVRPLIB

include("instanceCVRP.jl")
include("model.jl")


#global sumWeights = 0.0
global inst = createInstance("E-n76-k7") # E-n76-k7 E-n13-k4 E-n33-k4 

#print("sumWeights: ", sumWeights, "\n\n")
#exit(-1)

global vetReducedCosts = zeros(inst.n*inst.n)
global vetRoutesOutput = zeros(Cint, getSizeVetRoutesOutput())
global vetReducedCostsRoutes = zeros(getMaxNumberOfRoutes())
createModel(inst)

