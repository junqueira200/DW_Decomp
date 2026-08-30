include("packing.jl")

function build_model(data::DataArcVRP, app)

   A = arcs(data) 
   V = [i for i in 1:n(data)] 
   Q = veh_capacity(data)
   Vol_veic = veh_volume(data)

   useVolume  = true

   println("Q: ", Q)

   # Formulation
   vrptw = VrpModel()
   @variable(vrptw.formulation, x[a in A], Int)
   @objective(vrptw.formulation, Min, sum(c(data,a) * x[a] for a in A))
   @constraint(vrptw.formulation, indeg[i in V], sum(x[a] for a in A if a[2] == i) == 1.0)

   #println(vrptw.formulation)

   # Build the model directed graph G=(V1,A1)
   function linearIndex(i, j, tam)
           return i*tam + j
   end
   function buildgraph()

      v_source = v_sink = 0
      V1 = [i for i in 0:n(data)]

      L, U = lowerBoundNbVehicles(data), upperBoundNbVehicles(data) # multiplicity
      print("lowerBoundNbVehicles: ", lowerBoundNbVehicles(data), "\n\n")
      G = VrpGraph(vrptw, V1, v_source, v_sink, (L, U))

      epValues = [2E-5, 2E-7, 2E-9]
      dualVolumeTotal = getDualVolumeTotal()
      println("dualVolumeTotal: ", dualVolumeTotal)
      epValuesSize = length(epValues)
      numOfDualFunctions = epValuesSize * epValuesSize
      dictDualFeasiFunc = Dict{Int, Dict{Int, Float64}}()
      dictIndexRes_id =  Dict{Int, Int}()

      i = 0
      j = 0

      for ep0 in epValues

         j = 0
         for ep2 in epValues            
            setDualFeasibleFunction(ep0, ep2)

            dict = Dict{Int, Float64}()            
            index = linearIndex(i, j, length(epValues))

            dict[0] = 0.0
            
            for t in 0:n(data)
               dict[t] = getDualVolume(t)               
            end
         
            dictDualFeasiFunc[index] = dict
            #dictIndexRes_id[index] = add_resource!(G, main=false)

            j += 1
         end

         i += 1
      end

#      setDualFeasibleFunction(1.0/100.0, 1.0/150.0)

      #if app["enable_cap_res"]
      cap_res_id = add_resource!(G, main = true)
   
      if useVolume
         vol_res_id = add_resource!(G, main=false)
         #vol_dual_id = add_resource!(G, main=false)
      end
      #end
      

      for v in V1
         #if app["enable_cap_res"]
         set_resource_bounds!(G, v, cap_res_id, 0, Q)
         if useVolume
            set_resource_bounds!(G, v, vol_res_id, 0, Vol_veic)
            #set_resource_bounds!(G, v, vol_dual_id, 0.0, 1.0)

            for i_ in 0:(length(epValues)-1)
               for j_ in 0:(length(epValues)-1)
                  index = linearIndex(i_, j_, length(epValues))
                  #set_resource_bounds!(G, v, dictIndexRes_id[index], 0.0, dualVolumeTotal)
               
               end
            end
         end
         #end
         
      end

      for (i,j) in A
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])         
         set_arc_consumption!(G, arc_id, cap_res_id, d(data, j))

         if useVolume
            set_arc_consumption!(G, arc_id, vol_res_id, vol(data, j))

            for i_ in 0:(length(epValues)-1)
               for j_ in 0:(length(epValues)-1)
                  index = linearIndex(i_, j_, length(epValues))
                  #set_arc_consumption!(G, arc_id, dictIndexRes_id[index], dictDualFeasiFunc[index][j])
               
               end
            end

            #set_arc_consumption!(G, arc_id, vol_dual_id, getDualVolume(j))
         end
      end

      return G
   end

   G = buildgraph()
   add_graph!(vrptw, G)
   #println(G)

   set_vertex_packing_sets!(vrptw, [[(G,i)] for i in V])

   define_elementarity_sets_distance_matrix!(vrptw, G, [[c(data, (i, j)) for j in V] for i in V])

   add_capacity_cut_separator!(vrptw, [ ( [(G,i)], Float64(d(data, i)) ) for i in V], Float64(Q))

   set_branching_priority!(vrptw, "x", 1)

   return (vrptw, x)
end
