__precompile__(false)
module VRPTWSolverDemo
using VrpSolver, JuMP, ArgParse


include("data.jl")
include("model.jl")
include("solution.jl")
include("packing.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage="##### VRPSolver #####\n\n" *
                               "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
    @add_arg_table s begin
        "instance"
            help = "Instance file path"
        "--update"
            help = "Update the VrpSolver package"
            action = :store_true
        "--cfg", "-c"
            help = "Configuration file path"
            default = "$appfolder/../config/VRPTW_set_1.cfg"
        "--ub", "-u"
            help = "Upper bound (primal bound)"
            arg_type = Float64
            default = 10000000.0
        "--enable_cap_res", "-e"
            help = "Enable capacity resource"
            action = :store_true
        "--sol", "-s"
            help = "Solution file path (e.g., see sol/C203.sol)"
        "--out", "-o"
            help = "Path to write the solution found"
        "--tikz", "-t"
            help = "Path to write the TikZ figure of the solution found."
        "--nosolve", "-n"
            help = "Does not call the VRPSolver. Only to draw a given solution."
            action = :store_true
        "--batch", "-b"
            help = "batch file path"
        "--Oroloc3D"
            help = "Indicates the instance type is Oroloc3D"
            arg_type = Bool
            default  = true
        "--3dBasicPacking"
            help = "Indicates the instance type is 3dBasicPacking"
            arg_type = Bool
            default  = false

    end
    return parse_args(args_array, s)
end

function run_vrptw(app::Dict{String,Any})
    println("Application parameters:")
    for (arg, val) in app
        println("  $arg  =>  $(repr(val))")
    end
    flush(stdout)

    instance_name = split(basename(app["instance"]), ".")[1]

    data = createInstance(app["instance"], app["Oroloc3D"])

    if app["sol"] != nothing
        sol = readsolution(app["sol"])
        app["ub"] = (sol.cost < app["ub"]) ? sol.cost : app["ub"] # update the upper bound if necessary
    end

    solution_found = false
    if !app["nosolve"]
        (model, x) = build_model(data, app)
        optimizer = VrpOptimizer(model, app["cfg"], instance_name)
        set_cutoff!(optimizer, app["ub"])

        E = data.G′.A
        function mycallback()

            println("From mycallback")
            integer = true
            for (i,j) in E
                e = (i,j)
                value = get_value(optimizer, x[e])
                if abs(value - round(value)) > 1e-5
                    integer = false
                    break
                end
            end

            if integer
                println("Soluton is integer!")
                sol = getsolution(data, x, get_objective_value(optimizer), optimizer)
                for route in sol.routes
                    #println(route)
                    vet = Vector{Int32}()
                    push!(vet, 0)

                    for i in route
                        push!(vet, i)
                    end
                    
                    push!(vet, 0)
                    result = testRoute(vet)
                    if(result == 0)
                        for i in 1:(length(vet)-1)
                            arc = (vet[i], vet[i+1])
                            println("\t", arc)
                        end
                        println("Cut route")
                        
                        exit(-1)
                    end
                    #println(vet)

                end

                
                #add_dynamic_constr!(model.optimizer, [x[e]], [1.0], <=, 1.0, "edge_ub")
            end

        end 

        add_cut_callback!(model, mycallback, "mycallback")


        (status, solution_found) = optimize!(optimizer)
        if solution_found
            sol = getsolution(data, x, get_objective_value(optimizer), optimizer)
        end
    end

    println("########################################################")
    if solution_found || app["sol"] != nothing # Is there a solution?
        print_routes(sol)
        println("Cost $(sol.cost)")
        if app["out"] != nothing
            writesolution(app["out"], sol)
        end
        if app["tikz"] != nothing
            drawsolution(app["tikz"], data, sol) # write tikz figure
        end
    elseif !app["nosolve"]
        if status == :Optimal
            println("Problem infeasible")
        else
            println("Solution not found")
        end
    end
    println("########################################################")
end

function main(args)
    appfolder = dirname(@__FILE__)
    app = parse_commandline(args, appfolder)

    println("Oroloc3D: ", app["Oroloc3D"])
    #println("3dBasicPacking: ", app["3dBasicPacking"])

    isnothing(app) && return
    if app["batch"] != nothing
        for line in readlines(app["batch"])
            if isempty(strip(line)) || strip(line)[1] == '#'
                continue
            end
            args_array = [String(s) for s in split(line)]
            app_line = parse_commandline(args_array, appfolder)
            run_vrptw(app_line)
        end
    else
        run_vrptw(app)
    end
end



export main

end