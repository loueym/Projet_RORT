# MILP resolution 

using JuMP
using CPLEX

include("utils.jl")
include("data/taxe_grille_7x11.txt")

function get_milp_model(instance::Instance)
    # Model
    model = Model(CPLEX.Optimizer)
    
    # Variables
    @variable(model, T[a in eachindex(instance.A1)] >= 0)
    @variable(model, x1[a in eachindex(instance.A1), 1:instance.K] >= 0)
    @variable(model, x2[a in eachindex(instance.A2), 1:instance.K] >= 0)
    @variable(model, y[a in eachindex(instance.A1), k = 1:instance.K, l = 0:instance.demands[k]], Bin)
    @variable(model, z[a in eachindex(instance.A1), k = 1:instance.K, l = 0:instance.demands[k]] >= 0)
    @variable(model, lambda[1:instance.n, 1:instance.K])
    
    # Primal Constraints
    @constraint(model, [k = 1:instance.K, i in 1:instance.n; i != instance.origins[k] && i != instance.destinations[k]], 
                       sum(x1[a,k] for a in get_delta_plus_A1(instance, i)) + sum(x2[a,k] for a in get_delta_plus_A2(instance, i)) 
                       == sum(x1[a,k] for a in get_delta_minus_A1(instance, i)) + sum(x2[a,k] for a in get_delta_minus_A2(instance, i)))
    @constraint(model, [k in 1:instance.K], 
                       sum(x1[a,k] for a in get_delta_minus_A1(instance, instance.destinations[k])) 
                       + sum(x2[a,k] for a in get_delta_minus_A2(instance, instance.destinations[k])) 
                       == instance.demands[k])
    @constraint(model, [k in 1:instance.K], 
                       sum(x1[a,k] for a in get_delta_plus_A1(instance, instance.origins[k])) 
                       + sum(x2[a,k] for a in get_delta_plus_A2(instance, instance.origins[k])) 
                       == instance.demands[k])
    # Dual Constraints
    @constraint(model, [(i,j) in eachindex(instance.A1), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A1[i,j] + T[(i,j)])
    @constraint(model, [(i,j) in eachindex(instance.A2), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A2[i,j])
    @constraint(model, [k in 1:instance.K], lambda[instance.origins[k],k] == 0)
    # Primal-Dual Constraint
    @constraint(model, sum(l * z[a,k,l] for a in eachindex(instance.A1) for k in 1:instance.K for l in 0:instance.demands[k])
                       + sum(instance.A1[i,j] * x1[(i,j),k] for (i,j) in eachindex(instance.A1) for k in 1:instance.K)
                       + sum(instance.A2[i,j] * x2[(i,j),k] for (i,j) in eachindex(instance.A2) for k in 1:instance.K) 
                       <= sum(instance.demands[k] * lambda[instance.destinations[k], k] for k in 1:instance.K))
    # Linearization Constraints
    max_T_max = maximum(instance.T_max)
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K], x1[a,k] == sum(l * y[a,k,l] for l in 0:instance.demands[k]))
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K], sum(y[a,k,l] for l in 0:instance.demands[k]) == 1)
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K, l = 0:instance.demands[k]], z[a,k,l] <= y[a,k,l] * max_T_max)
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K, l = 0:instance.demands[k]], z[a,k,l] <= T[a])
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K, l = 0:instance.demands[k]], z[a,k,l] >= T[a] - (1 - y[a,k,l]) * max_T_max)
    

    # Objective
    @objective(model, Max, sum(l * z[a,k,l] for a in eachindex(instance.A1) for k = 1:instance.K for l in 0:instance.demands[k]))
    return model
end

function get_compact_model(instance::Instance)
    # Model
    model = Model(CPLEX.Optimizer)
    
    # Variables
    @variable(model, T[a in eachindex(instance.A1)] >= 0)
    @variable(model, y1[a in eachindex(instance.A1), k = 1:instance.K], Bin)
    @variable(model, y2[a in eachindex(instance.A2), k = 1:instance.K], Bin)
    @variable(model, z[a in eachindex(instance.A1), k = 1:instance.K] >= 0)
    @variable(model, lambda[1:instance.n, 1:instance.K])
    
    # Primal Constraints
    @constraint(model, [k = 1:instance.K, i in 1:instance.n; i != instance.origins[k] && i != instance.destinations[k]], 
                       sum(y1[a,k] for a in get_delta_plus_A1(instance, i)) + sum(y2[a,k] for a in get_delta_plus_A2(instance, i)) 
                       == sum(y1[a,k] for a in get_delta_minus_A1(instance, i)) + sum(y2[a,k] for a in get_delta_minus_A2(instance, i)))
    @constraint(model, [k in 1:instance.K], 
                       sum(y1[a,k] for a in get_delta_minus_A1(instance, instance.destinations[k])) 
                       + sum(y2[a,k] for a in get_delta_minus_A2(instance, instance.destinations[k])) 
                       == 1) 
    @constraint(model, [k in 1:instance.K], 
                       sum(y1[a,k] for a in get_delta_plus_A1(instance, instance.origins[k])) 
                       + sum(y2[a,k] for a in get_delta_plus_A2(instance, instance.origins[k])) 
                       == 1)
    # Dual Constraints
    @constraint(model, [(i,j) in eachindex(instance.A1), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A1[i,j] + T[(i,j)])
    @constraint(model, [(i,j) in eachindex(instance.A2), k in 1:instance.K], lambda[j,k] - lambda[i,k] <= instance.A2[i,j])
    @constraint(model, [k in 1:instance.K], lambda[instance.origins[k],k] == 0)
    # Primal-Dual Constraint
    @constraint(model, [k in 1:instance.K], 
                       sum(instance.A1[i,j] * y1[(i,j),k] + z[(i,j),k] for (i,j) in eachindex(instance.A1))
                       + sum(instance.A2[i,j] * y2[(i,j),k] for (i,j) in eachindex(instance.A2)) 
                       == lambda[instance.destinations[k], k])
    # Linearization Constraints
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K], z[a,k] <= y1[a,k] * instance.T_max[k])
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K], z[a,k] <= T[a])
    max_T_max = maximum(instance.T_max)
    @constraint(model, [a in eachindex(instance.A1), k = 1:instance.K], z[a,k] >= T[a] - (1 - y1[a,k]) * max_T_max)

    # Objective
    @objective(model, Max, sum(instance.demands[k] * z[a,k] for a in eachindex(instance.A1) for k = 1:instance.K))

    # write_to_file(model, "test_model.lp")
    return model
end

function solve_milp_model!(instance::Instance, model::Model, compact::Bool=true)
    set_time_limit_sec(model, 60)
    optimize!(model)
    feasible_solution_found = primal_status(model) == MOI.FEASIBLE_POINT
    is_optimal = termination_status(model) == MOI.OPTIMAL

    if feasible_solution_found || is_optimal
        vTaxes = JuMP.value.(model[:T])
        # display(vTaxes)
        # println()
        # lambdas = JuMP.value.(model[:lambda])
        # for k in 1:instance.K
        #     cost_k = round(lambdas[instance.destinations[k], k])
        #     println("Cout du chemin pour commodité $k : $(cost_k)")
        # end
        # println()
        # if compact
        #     y1 = JuMP.value.(model[:y1])
        #     y2 = JuMP.value.(model[:y2])
        #     for k in 1:instance.K
        #         for a in eachindex(instance.A1)
        #             if y1[a,k] > 1- 1e-5
        #                 println("Commodité $k sur arc $a de cout $(instance.A1[a]) et de taxe $(round(vTaxes[a]))")
        #             end
        #         end
        #         for a in eachindex(instance.A2)
        #             if y2[a,k] > 1- 1e-5
        #                 println("Commodité $k sur arc $a de cout $(instance.A2[a])")
        #             end
        #         end
        #         println()
        #     end
        # else
        #     x1 = JuMP.value.(model[:x1])
        #     x2 = JuMP.value.(model[:x2])
        #     for k in 1:instance.K
        #         for a in eachindex(instance.A1)
        #             if x1[a,k] > 1- 1e-5
        #                 println("Commodité $k sur arc $a de cout $(instance.A1[a]) et de taxe $(round(vTaxes[a]))")
        #             end
        #         end
        #         for a in eachindex(instance.A2)
        #             if x2[a,k] > 1- 1e-5
        #                 println("Commodité $k sur arc $a de cout $(instance.A2[a])")
        #             end
        #         end
        #         println()
        #     end
        # end
        obj = JuMP.objective_value(model)
        println("Profit from taxes : $(obj)")
        z_star = JuMP.value.(model[:z])

        return vTaxes, obj, z_star
    end
end


function print_shortest_paths(z_star, instance::Instance)
    # A, K = size(z_star)
    # println(A, K)
    println(z_star)
    # for a in 1:A
    #     if sum(z_star[a, :]) > 0
    #         print("On arc ")
    #     for k in 1:K

    #     end
    # end
end

# Test
test_instance = Instance(K, n_k, n, A_1, A_2)
println(test_instance)
# test_model = get_milp_model(test_instance)
# show(test_model)
# println()
# test_solution1, test_obj1, test_z = solve_milp_model!(test_instance, test_model, false)
# println()
# println(test_solution)

test_model = get_compact_model(test_instance)
show(test_model)
println()
test_solution2, test_obj2, test_z = solve_milp_model!(test_instance, test_model)
println()
println(test_solution2)
println()
# print_shortest_paths(test_z, test_instance)

# test_model = get_milp_model(test_instance)
# show(test_model)
# test_solution = solve_milp_model!(test_model)
# println(test_solution)
