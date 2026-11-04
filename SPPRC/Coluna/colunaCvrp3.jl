using JuMP, BlockDecomposition, Coluna, Gurobi
const BD  = BlockDecomposition
const MOI = JuMP.MOI
using CVRPLIB

include("instanceCVRP.jl")
include("model.jl")


#global sumWeights = 0.0
global inst = createInstance("E-n76-k7") # E-n76-k7 E-n13-k4 E-n33-k4 

print("sumWeights: ", sumWeights, "\n\n")
#exit(-1)

global vetReducedCosts = zeros(inst.n*inst.n)
global vetRoutesOutput = zeros(Cint, getSizeVetRoutesOutput())
global vetReducedCostsRoutes = zeros(getMaxNumberOfRoutes())
createModel(inst)

exit(-1)


# ---------------- Dados da instância (1 = depósito) ----------------
coord  = [(0,0), (2,4), (5,1), (6,6), (1,8), (7,3)]
demand = [0, 4, 3, 5, 4, 2]
Q, K   = 10, 3                        # capacidade e nº máximo de veículos

V         = 1:length(coord)
customers = 2:length(coord)
arcs      = [(i, j) for i in V for j in V if i != j]   # arcos dirigidos
cost      = Dict(a => hypot(coord[a[1]][1]-coord[a[2]][1],
                            coord[a[1]][2]-coord[a[2]][2]) for a in arcs)

# ---------------- Enumeração prévia das rotas viáveis ----------------
# Rota = depósito -> sequência de clientes distintos (carga <= Q) -> depósito
function enumerate_routes(customers, demand, Q)
    routes = Vector{Vector{Int}}()
    function extend(path, load)
        isempty(path) || push!(routes, [1; path; 1])
        for c in customers
            c in path && continue
            load + demand[c] <= Q || continue
            extend([path; c], load + demand[c])
        end
    end
    extend(Int[], 0)
    return routes
end

routes     = enumerate_routes(customers, demand, Q)
route_arcs = [[(r[k], r[k+1]) for k in 1:length(r)-1] for r in routes]
@info "Rotas viáveis enumeradas" n = length(routes)

# ---------------- Configuração do Coluna ----------------

coluna = optimizer_with_attributes(
    Coluna.Optimizer,
    "params" => Coluna.Params(
        solver = Coluna.Algorithm.TreeSearchAlgorithm(
            conqueralg = Coluna.Algorithm.ColCutGenConquer(
                colgen = Coluna.Algorithm.ColumnGeneration(
                    smoothing_stabilization = 0.7   # em [0, 1]
                )
            )
        )
    ),
    "default_optimizer" => GLPK.Optimizer,
)
model = BlockModel(coluna)

@axis(VehicleTypes, [1])

@variable(model, x[v in VehicleTypes, a in arcs], Bin)

# Mestre: cada cliente é visitado exatamente uma vez
@constraint(model, visit[i in customers],
    sum(x[1, a] for a in arcs if a[2] == i) == 1)

@objective(model, Min, sum(cost[a] * x[1, a] for a in arcs))

# ---------------- Callback do subproblema ----------------
function pricing_callback(cbdata)
    v = BD.callback_spid(cbdata, model)

    # 1) Recupera o custo reduzido de cada arco
    rc = Dict(a => BD.callback_reduced_cost(cbdata, x[v, a]) for a in arcs)

    # 2) Custo reduzido de cada rota enumerada = soma dos seus arcos;
    #    guarda a menor
    best_k, best_val = 0, Inf
    for (k, ra) in enumerate(route_arcs)
        val = sum(rc[a] for a in ra)
        if val < best_val
            best_k, best_val = k, val
        end
    end

    # 3) Submete a melhor rota como coluna e o dual bound (pricing exato)
    best = route_arcs[best_k]
    MOI.submit(model, BD.PricingSolution(cbdata), best_val,
               [x[v, a] for a in best], ones(length(best)))
    MOI.submit(model, BD.PricingDualBound(cbdata), best_val)
end

# ---------------- Decomposição de Dantzig-Wolfe ----------------
@dantzig_wolfe_decomposition(model, decomposition, VehicleTypes)
subproblems = BD.getsubproblems(decomposition)
BD.specify!.(subproblems,
             lower_multiplicity = 0,
             upper_multiplicity = K,
             solver = pricing_callback)

optimize!(model)

# ---------------- Recuperação da solução ----------------
@info "Resultado" termination_status(model) primal_status(model)

if primal_status(model) == MOI.FEASIBLE_POINT
    tol = 1e-6

    # Arcos ativos (a solução final do branch-and-price é inteira)
    active = [a for a in arcs if value(x[1, a]) > 1 - tol]
    frac   = [a for a in arcs if tol < value(x[1, a]) < 1 - tol]
    isempty(frac) || @warn "Arcos fracionários na solução final" frac

    # Sucessor de cada cliente e primeiro cliente de cada rota
    succ   = Dict{Int,Int}()
    starts = Int[]
    for (i, j) in active
        i == 1 ? push!(starts, j) : (succ[i] = j)
    end

    # Reconstrói cada rota: depósito -> clientes -> depósito
    solution_routes = Vector{Vector{Int}}()
    for s in starts
        route = [1, s]
        while route[end] != 1
            haskey(succ, route[end]) ||
                error("Arco de saída ausente no nó $(route[end])")
            push!(route, succ[route[end]])
        end
        push!(solution_routes, route)
    end

    # Custo de cada rota (sem acumulador global)
    route_costs = [sum(cost[(r[t], r[t+1])] for t in 1:length(r)-1)
                   for r in solution_routes]

    # Impressão e verificação
    println("\n========== Solução ==========")
    for (k, r) in enumerate(solution_routes)
        load = sum(demand[c] for c in r)
        id   = findfirst(==(r), routes)   # índice na lista enumerada
        println("Rota $k: ", join(r, " → "),
                " | carga = $load/$Q",
                " | custo = $(round(route_costs[k]; digits = 3))",
                " | rota enumerada #$id")
    end
    println("Custo total das rotas:   ", round(sum(route_costs); digits = 4))
    println("Valor objetivo (Coluna): ", round(objective_value(model); digits = 4))

    # Verificações de sanidade
    visited = sort(vcat([r[2:end-1] for r in solution_routes]...))
    @assert visited == collect(customers) "Nem todos os clientes foram visitados exatamente uma vez"
    @assert all(sum(demand[c] for c in r) <= Q for r in solution_routes) "Capacidade violada"
    @assert length(solution_routes) <= K "Excede o número de veículos"
    @assert isapprox(sum(route_costs), objective_value(model); atol = 1e-4) "Custo reconstruído difere do objetivo"
end
