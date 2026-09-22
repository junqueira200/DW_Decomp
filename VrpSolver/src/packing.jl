include("data.jl")

function testRouteCapVol(route::Vector{Int32})
    return ccall((:testRouteCapVol, LIB),
                 Cint,
                 (Ptr{Cint}, Cint),
                 route, 
                 length(route))
end

function getDistanceRoute(route::Vector{Int32})
    return ccall((:getDistanceRoute, LIB),
                Cdouble,
                (Ptr{Cint}, Cint),
                route, 
                length(route))
end

function testRoute(route::Vector{Int32})
    
    println("Route: ", route)

    result = ccall(
                (:testRoute, LIB),
                Cint,
                (Ptr{Cint}, Cint, Cint, Cint),
                route, 
                length(route),
                0,
                0                
                )

    println(route, ": ", result)
    return result

end


function setDualFeasibleFunction(ep0::Float64, ep2::Float64)
    ccall(
        (:setDualFeasibleFunction, LIB),
        Cvoid,
        (Cdouble, Cdouble),
        ep0,
        ep2
        ) 
end

function getDualVolume(cust::Int)
    return ccall(
        (:getDualVolume, LIB),
        Cdouble,
        (Cint,),
        cust
    )
end

function getDualVolumeTotal()
    return ccall(
        (:getDualVolumeTotal, LIB),
        Cdouble,
        ()
    )
end

function omp_get_wtime()
    return ccall(
        (:omp_get_wtime, LIB),
        Cdouble,
        ()
    )
end

function getUseExactPacking()
    return ccall((:getUseExactPacking, get(ENV, "BAPCOD_RCSP_LIB", "")), Cint, ())
    
end

function setFalseHeuristicWorks()
    println("setFalseHeuristicWorks")
    ccall((:getPtrRoutes, get(ENV, "BAPCOD_RCSP_LIB", "")), Cvoid, ())
end


function getMasterRutesIsNull()
    len = ccall((:getSizePtrRoutes, get(ENV, "BAPCOD_RCSP_LIB", "")), Cint, ())

    return len == 0
end
function getMasterRutes()    
    
    len = ccall((:getSizePtrRoutes, get(ENV, "BAPCOD_RCSP_LIB", "")), Cint, ())
    ptr = ccall((:getPtrRoutes, get(ENV, "BAPCOD_RCSP_LIB", "")), Ptr{Cint}, ())

    return unsafe_wrap(Array, ptr, len, own=false)

end

function routeIsInFeasibleSet(route::Vector{Int32})
    
    println("Route: ", route)

    result = ccall(
                (:routeIsInFeasibleSet, LIB),
                Cint,
                (Ptr{Cint}, Cint),
                route, 
                length(route), 
                )

    
    if result >= 1
        println(route, " is feasible")
    end
    return result

end

function routeIsInNotfeasibleSet(route::Vector{Int32})
        
    result = ccall(
                (:routeIsInNotfeasibleSet, LIB),
                Cint,
                (Ptr{Cint}, Cint),
                route, 
                length(route), 
                )

    if result >= 1
        println(route, " is infeasible")
    end

    return result

end


function heuristicPacking(route::Vector{Int32})
    
    result = ccall(
                (:heuristicPacking, LIB),
                Cint,
                (Ptr{Cint}, Cint),
                route, 
                length(route), 
                )

    if result >= 0
        println("heuristicPacking: ", route)

    else
        print("heuristicPacking: ", route, " faild")
    end

    return result
    
end

function exactPacking(route::Vector{Int32})
    
    result = ccall(
                (:exactPacking, LIB),
                Cint,
                (Ptr{Cint}, Cint),
                route, 
                length(route), 
                )

    if result >= 1
        println("exactPacking: ", route)

    else
        print("exactPacking: ", route, " faild")
    end

    return result
    
end

