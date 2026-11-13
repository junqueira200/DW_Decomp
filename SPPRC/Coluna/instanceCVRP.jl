using CVRPLIB

const LIB = "../build/libSPPRC.dylib"

function setInstance(n)
    ccall((:setInstance, LIB),
           Cvoid,
           (Cint,),
            n)
end

function createArcs(n, weights)

    arc_weights = Dict{Tuple{Int, Int}, Int}()
    arcs = Vector{Tuple{Int, Int}}()
    global sumWeights = 0.0

    for i in 0:n-1
        for j in 0:n-1
        
            arc_weights[(i,j)] = weights[i+1,j+1]
            sumWeights += weights[i+1,j+1]

            if i != j
                push!(arcs, (i, j))
            end
            
        end
    end

    sumWeights = sumWeights*2

    return arc_weights, arcs

end

function createDemands(n, demands)
    
    dictDemand = Dict{Int, Int}()

    for i in 1:n
        dictDemand[i-1] = demands[i]
    end
    return dictDemand
end

struct CvrpData

    n::Int
    arc_weights::Dict{Tuple{Int, Int}, Int}
    arcs::Vector{Tuple{Int, Int}}
    demand::Dict{Int, Int}
    customers::UnitRange{Int64}
    numberOfVehicles::Int
    capacity::Int

end

function setDistance(inst::CvrpData)
    
    for tuple in inst.arc_weights
        arc    = tuple[1]
        weight = tuple[2]

        if arc[1] != arc[2]
            #println(arc, ": ", weight)
            ccall((:setDistance, LIB),
                    Cvoid,
                    (Cint, Cint, Cint,),
                    arc[1],
                    arc[2],
                    weight
                 )
        end
    end

end

function setDemand(inst::CvrpData)
    
    for tuple in inst.demand
        ccall((:setDemand, LIB),
                Cvoid,
                (Cint, Cint,),
                tuple[1],
                tuple[2])
    end

end

function setCapacityVehicle(inst::CvrpData)
    ccall((:setCapacityVehicle, LIB), Cvoid, (Cint,), inst.capacity)
end

function setVrpLabelingSubProb()

    ccall((:setVrpLabelingSubProb, LIB), Cvoid, ())
            
end

function setNumberOfVehicles(k::Int)
    ccall((:setNumberOfVehicles, LIB), Cvoid, (Cint,), k)
end

function getSizeVetRoutesOutput()
    return ccall((:getSizeVetRoutesOutput, LIB), Cint, ())
end

function getMaxNumberOfRoutes()
    return ccall((:getMaxNumberOfRoutes, LIB), Cint, ())
end

function createInstance(instanceStr::String)
    cvrp, _, _ = readCVRPLIB(instanceStr)
    n = maximum(cvrp.customers) - 1

    arc_weights, arcs = createArcs(n, cvrp.weights)
    dictDemand = createDemands(n, cvrp.demand)
    sumDemands = Float64(sum(cvrp.demand))
    k = Int(round(sumDemands/Float64(cvrp.capacity), RoundUp))
    println("numberOfVehicles: ", k)


    inst = CvrpData(n, arc_weights, arcs, dictDemand, 1:n-1, k, cvrp.capacity)
    #println(cvrp)

    setInstance(n)

    setDistance(inst)
    setDemand(inst)
    setCapacityVehicle(inst)
    setNumberOfVehicles(inst.numberOfVehicles)

    setVrpLabelingSubProb()
    println("After")
    return inst
end