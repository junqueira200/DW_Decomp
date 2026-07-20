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