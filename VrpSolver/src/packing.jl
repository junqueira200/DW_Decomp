include("data.jl")

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