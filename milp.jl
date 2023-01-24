# MILP resolution 

using JuMP
using CPLEX

include("instance.jl")
include("data/taxe_grille_2x3.txt")

function get_milp_model(instance::Instance)
    # Model
    model = Model(CPLEX.Optimizer)
    
    # Variables
    @variable(model, T[(i,j) in eachindex(instance.A1)] >= 0)
    @variable(model, x1[(i,j) in eachindex(instance.A1), 1:instance.K] >= 0)
    @variable(model, x2[(i,j) in eachindex(instance.A2), 1:instance.K] >= 0)
    @variable(model, y[(i,j) in eachindex(instance.A1), k = 1:instance.K, l = 1:instance.demands[k]], Bin)
    @variable(model, z[(i,j) in eachindex(instance.A1), k = 1:instance.K, l = 1:instance.demands[k]] >= 0)
    @variable(model, lambda[1:instance.n, 1:instance.K])
    
    # Primal Constraints
    @constraint(model, [i in 1:instance.n], sum(x[i,k] for k in 1:instance.K) == 1) # every point in a cluster
    @constraint(model, [k in 1:instance.K], sum(instance.w_v[i] * x[i,k] for i in 1:instance.n) <= instance.B) # cluster weights
    @constraint(model, [i = 1:instance.n, j = (i+1):instance.n, k in 1:instance.K], 1 + y[i,j] >= x[i,k] + x[j,k]) # linearization constraint
    # Dual Constraints
    @constraint(model, [(i,j) in eachindex(instance.A1), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A1[i,j] + T[i,j])
    @constraint(model, [(i,j) in eachindex(instance.A2), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A2[i,j])
    @constraint(model, [k in 1:instance.K], lambda[instance.origins[k],k] == 0)
    # Primal-Dual Constraint
    @constraint(model, sum() + sum() <= sum(instance.demands[k] * lambda[instance.destinations[k], k] for k in 1:instance.K))
    # Linearization Constraints
    @constraint(model, [(i,j) in eachindex(instance.A1), k = 1:instance.K], x1[i,j,k] == sum(l * y[i,j,k,l] for l in 1:instance.demands[k]))
    @constraint(model, [(i,j) in eachindex(instance.A1), k = 1:instance.K, l = 1:instance.demands[k]], z[i,j,k,l] <= y[i,j,k,l] * instance.T_max)
    @constraint(model, [(i,j) in eachindex(instance.A1), k = 1:instance.K, l = 1:instance.demands[k]], z[i,j,k,l] <= T[i,j])
    @constraint(model, [(i,j) in eachindex(instance.A1), k = 1:instance.K, l = 1:instance.demands[k]], z[i,j,k,l] >= T[i,j] - (1 - y[i,j,k,l]) * instance.T_max)
    

    # Objective
    @objective(model, Min, sum(instance.distances[i,j] * y[i,j] for i in 1:instance.n, j in (i+1):instance.n))
    return model
end

# Test
test_instance = Instance(K, n_k, n, A_1, A_2)