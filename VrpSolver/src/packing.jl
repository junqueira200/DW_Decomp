include("data.jl")

function testRoute(route::Vector{Int32})
    
    println("Route: ", route)

    result = ccall(
                (:testRoute, LIB),
                Cint,
                (Ptr{Cint}, Cint),
                route, 
                length(route)                
                )

    println(route, ": ", result)
    return result

end