const LIB = "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/build/2L_SDVRP/lib2L_SDVRP.dylib"

#const inst = "/Users/igor/Documents/[08] Projetos/[16]_DW_Decomp/2L_SDVRP/Instances/Oroloc3D_2/instances/S10_240T.txt"

vet = [0, 0]
cVet = Cint.(vet)

#println(cVet)

mutable struct Vertex
   id_vertex::Int
   #pos_x::Float64
   #pos_y::Float64
   #service_time::Int
   demand::Int
   volume::Float64
   #start_tw::Int
   #end_tw::Int

end

# Directed graph
mutable struct InputGraph
   V′::Vector{Vertex} # set of vertices
   A::Vector{Tuple{Int64,Int64}} # set of edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each arc
   #time::Dict{Tuple{Int64,Int64},Float64} # time for each arc
end

mutable struct DataArcVRP
   n::Int
   G′::InputGraph
   Q::Float64 # vehicle capacity
   K::Int #num vehicles available
   Vol::Float64  # vehicle volume
end


function ini_3D_Packing(instt, oroloc3D)

    ccall(
        (:ini_3D_Packing, LIB),
        Cvoid,
        (Cstring, Cint),
        instt,
        oroloc3D
    )
    
end

function omp_get_wtime()
    return ccall(
        (:omp_get_wtime, LIB),
        Cdouble,
        ()
    )
end

function getNumberOfCustoms()
    return ccall(
                (:getNumberOfCustoms, LIB),
                Cint,
                (),                
            )
end

function getNumberOfTrucks()
    return ccall(
                (:getNumberOfTrucks, LIB),
                Cint,
                (),                
            )        
end

function getDemandFromCustomr(id::Int)
    return ccall(
                (:getDemandFromCustomr, LIB),
                Cint,
                (Cint,),
                id                
            )   
end

function getVolumeFromCustomr(id::Int)
    
    #println("id: ", id)

    return ccall(
                (:getVolumeFromCustomr, LIB),
                Cdouble,
                (Cint,),
                id                
            )   
end


function getDistance(i::Int, j::Int)
    # convert from mm to km
    return ccall(
                (:getDistance, LIB),
                Cdouble,
                (Cint,Cint,),
                i,
                j           
            )       
end

function getVehicleCapacity()
    return ccall(
                (:getVehicleCapacity, LIB),
                Cint,
                (),                
            )      
end

function getVehicleVolume()
    return ccall(
                (:getVehicleVolume, LIB),
                Cdouble,
                (),                
            )
end

function createInstance(instt, oroloc3D)
    
    ini_3D_Packing(instt, oroloc3D)


    numberOfCustoms = getNumberOfCustoms()
    numberOfTrucks  = getNumberOfTrucks()
    vehicleCapacity = getVehicleCapacity()
    vehicleVolume = getVehicleVolume()

    println("vehicleVolume: ", getVehicleVolume())

    demand = Dict{Int64,Int64}()
    volume = Dict{Int64,Float64}()
    dist   = Dict{Tuple{Int,Int}, Float64}()

    #println("numberOfTrucks: ", numberOfTrucks)

    sumDemand = 0
    sumVol    = 0.0

    for i in 0:numberOfCustoms-1
        demand[i] = getDemandFromCustomr(i)
        volume[i] = getVolumeFromCustomr(i)

        sumDemand += demand[i]
        sumVol    += volume[i]
    end
    
    minNumberOfTrucksDemand = Int(round(Float64(sumDemand)/Float64(vehicleCapacity)))
    minNumberOfTrucksVolume = Int(round(sumVol/Float64(vehicleVolume)))

    minNumberOfTrucksAux = max(minNumberOfTrucksDemand, minNumberOfTrucksVolume)
    #numberOfTrucks       = min(numberOfTrucks, minNumberOfTrucksAux)

    print("\nMin number of trucks (demand): ", minNumberOfTrucksDemand)
    print("Min number of trucks (volume): ", minNumberOfTrucksVolume, "\n\n")
    print("numberOfTrucks: ", numberOfTrucks, "\n\n")

    for i in 0:numberOfCustoms-1
        for j in 0:numberOfCustoms-1
            if i != j            
                dist[(i,j)] = getDistance(i, j)
            end
        end
    end

    println("volume: ", volume)

    println("demand: ", demand)
    println("*********************\n\n")
    println("vehicleCapacity: ", vehicleCapacity)
    #println("Dist: ", dist)
    #println("*********************\n\n")

    #println("vehicleCapacity: ", vehicleCapacity)
    #println("Demand: ", demand)

    vetVertex = Vector{Vertex}()
    vetEdges  = Tuple{Int64,Int64}[]

    for i in 0:numberOfCustoms-1
        push!(vetVertex, Vertex(i, demand[i], volume[i]))
        
        for j in 0:numberOfCustoms-1
            if i != j
                push!(vetEdges, (i,j))
            end
        end

    end



    #println(vetVertex)
    #println(vetEdges)
    DataArcVRP(numberOfCustoms-1, InputGraph(vetVertex, vetEdges, dist), vehicleCapacity, numberOfTrucks, vehicleVolume)

end

arcs(data::DataArcVRP) = data.G′.A # return set of arcs
function c(data,a) 
   if !(haskey(data.G′.cost, a)) 
      return Inf
   end
   return data.G′.cost[a] 
end

n(data::DataArcVRP) = data.n # return number of requests
d(data::DataArcVRP, i) = data.G′.V′[i+1].demand # return demand of i
veh_capacity(data::DataArcVRP) = Int(data.Q)
veh_volume(data::DataArcVRP) = data.Vol
vol(data::DataArcVRP, i) = data.G′.V′[i+1].volume

function lowerBoundNbVehicles(data::DataArcVRP) 
    sumVol = 0.0
    sumDemand = 0
    for it in data.G′.V′
        sumVol += it.volume
        sumDemand += it.demand
    end

    return Int(max(round(sumVol/veh_volume(data)), round(sumDemand/veh_capacity(data))))
end

function upperBoundNbVehicles(data::DataArcVRP) 
   return data.K
end


#data = createInstance(inst)
#print(data)


#println(inst)