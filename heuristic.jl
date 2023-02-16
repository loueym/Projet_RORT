# Heuristic resolution

using JuMP
using CPLEX
using Random

include("instance.jl")

function compute_shortest_path_with_pl(instance::Instance, journey_index::Int64, taxes)
    # to finish 

    # Model
    model = Model(CPLEX.Optimizer) 
    
    # Variables
    @variable(model, x1[(i,j) in eachindex(instance.A1)] >= 0)
    @variable(model, x2[(i,j) in eachindex(instance.A2)] >= 0)

    # add constraints

    # Objective
    @objective(
        model, 
        Min, 
        (
            sum(x1[i, j]*(instance.A1[i, j] + get(taxes, (i, j), 0)) for (i, j) in keys(instance.A1)) 
            + sum(x2[i, j]*instance.A2[i, j] for (i, j) in keys(instance.A2))
        )
    )

end


function compute_paths_through_e(instance::Instance, e::Tuple{Int, Int}, taxes::Array{Float64, 2})::Dict{}
    # returns a dictionnary with couples (origin, destination) as keys, and the max tax for (ori, des) on edge e as value
    inf = 100000
    adj_mat = build_adj_mat(instance.n, instance.A1, instance.A2, taxes)
    paths_through_e = Dict{}()
    for i in 1:instance.K
        ori = instance.origins[i]
        des = instance.destinations[i]
        dijk = dijkstra_shortest_paths(instance.g, ori, adj_mat)
        shortest_path = enumerate_paths(dijk, des)
        if e[1] in shortest_path && e[2] in shortest_path
            index1 = findfirst(item -> item == e[1], shortest_path)
            index2 = findfirst(item -> item == e[2], shortest_path)
            if index1 + 1 == index2
                # if the shortest_path visits e[2] right after e[1], we add it to our list
                dist = dijk.dists[des]
                cost, adj_mat[e[1], e[2]] = adj_mat[e[1], e[2]], inf
                next_dist = dijkstra_shortest_paths(instance.g, ori, adj_mat).dists[des]
                adj_mat[e[1], e[2]] = cost
                # the maximum tax is the gap between the cost of the second best path and the first one
                paths_through_e[(ori, des)] = next_dist - dist
            end
        end
    end
    return paths_through_e
end


function compute_obj_value(instance::Instance, taxes::Array{Float64, 2})
    val_before_taxes = 0
    val_after_taxes = 0
    n = instance.n
    adj_mat = build_adj_mat(n, instance.A1, instance.A2, zeros(n, n))
    adj_mat_taxes = build_adj_mat(n, instance.A1, instance.A2, taxes)
    for i in 1:instance.K
        ori = instance.origins[i]
        des = instance.destinations[i]
        demand = instance.demands[i]
        dist = dijkstra_shortest_paths(instance.g, ori, adj_mat).dists[des]
        dist_taxes = dijkstra_shortest_paths(instance.g, ori, adj_mat_taxes).dists[des]
        val_before_taxes += dist * demand
        val_after_taxes += dist_taxes * demand
    end
    return val_after_taxes - val_before_taxes
end


function initial_sol(instance::Instance, start_from_Tmax::Bool)
    n = instance.n
    taxes = zeros(n, n)
    if start_from_Tmax
        inf_mat(taxes, maximum(instance.T_max))
    end

    for e in shuffle!(collect(keys(instance.A1)))
        current_tax = 0
        paths_through_e = compute_paths_through_e(instance, e, taxes)
        possible_taxes = sort!(collect(values(paths_through_e)), rev=start_from_Tmax)
        previous_obj = 0 
        obj = compute_obj_value(instance, taxes)
        iter = 1
        while obj >= previous_obj && iter <= length(possible_taxes)
            # while we improve the objective, we raise the tax if start_from_Tmax is false,
            # and we lower it if start_from_Tmax=true
            previous_tax = current_tax
            current_tax = possible_taxes[iter]
            if current_tax != previous_tax
                previous_obj = obj
                taxes[e[1], e[2]] = current_tax
                obj = compute_obj_value(instance, taxes)
                if obj < previous_obj
                    # undo changes
                    taxes[e[1], e[2]] = previous_tax
                end
            end
            iter += 1
        end
    end
    return taxes
end

function heuristic(instance::Instance)
    taxes_from_zero = initial_sol(instance, false)
    value_from_zero = compute_obj_value(test_instance, taxes_from_zero)
    taxes_from_Tmax = initial_sol(instance, true)
    value_from_Tmax = compute_obj_value(test_instance, taxes_from_Tmax)
    if value_from_Tmax > value_from_zero
        return value_from_Tmax, taxes_from_Tmax
    end
    return value_from_zero, taxes_from_zero
end