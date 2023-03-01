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


function compute_paths_through_a(instance::Instance, arc::Tuple{Int, Int}, taxes::Array{Float64, 2})::Dict{}
    # returns a dictionnary with couples (origin, destination) as keys, and the max tax for (ori, des) on the arc as value
    n = instance.n
    n = instance.n
    adj_mat = build_dist_mat(instance.n, instance.A1, instance.A2, taxes)
    ad_mat_without_taxed_arcs = build_dist_mat(instance.n, instance.A1, instance.A2, INF .* ones(n, n))
    ad_mat_without_taxed_arcs = build_dist_mat(instance.n, instance.A1, instance.A2, INF .* ones(n, n))
    paths_through_a = Dict{}()
    for i in 1:instance.K
        ori = instance.origins[i]
        des = instance.destinations[i]
        dijk = dijkstra_shortest_paths(instance.g, ori, adj_mat)
        shortest_path = enumerate_paths(dijk, des)
        n1, n2 = arc
        if n1 in shortest_path && n2 in shortest_path
            index1 = findfirst(item -> item == n1, shortest_path)
            index2 = findfirst(item -> item == n2, shortest_path)
            if index1 + 1 == index2
                # if the shortest_path visits n2 right after n1, we add it to our list
                dist = dijk.dists[des]
                # cost, adj_mat[n1, n2] = adj_mat[n1, n2], inf
                next_dist = dijkstra_shortest_paths(instance.g, ori, ad_mat_without_taxed_arcs).dists[des]
                # adj_mat[n1, n2] = cost
                # cost, adj_mat[n1, n2] = adj_mat[n1, n2], inf
                next_dist = dijkstra_shortest_paths(instance.g, ori, ad_mat_without_taxed_arcs).dists[des]
                # adj_mat[n1, n2] = cost
                # the maximum tax is the gap between the cost of the second best path and the first one
                paths_through_a[(ori, des)] = next_dist - dist
            end
        end
    end
    return paths_through_a
end


function compute_obj_value(instance::Instance, taxes::Array{Float64, 2}, show_result::Bool)
    val = 0
    n = instance.n
    adj_mat = build_dist_mat(n, instance.A1, instance.A2, taxes)
    for i in 1:instance.K
        ori = instance.origins[i]
        des = instance.destinations[i]
        demand = instance.demands[i]
        dijk = dijkstra_shortest_paths(instance.g, ori, adj_mat)
        shortest_path = enumerate_paths(dijk, des)
        prev_n = shortest_path[1]
        if show_result
            println("Shortest path from ", ori, " to ", des, " is ", shortest_path, " with cost ", dijk.dists[des])
        end
        taxes_along_path = 0
        for n in shortest_path[2:length(shortest_path)]
            if show_result
                println("for arc (", prev_n, ", ", n, "), taxes = ", taxes[prev_n, n])
            end
            taxes_along_path += taxes[prev_n, n] * demand
            prev_n = n
        end
        val += taxes_along_path
        if show_result
            println("-> taxes along shortest path are ", taxes_along_path)
            println()
        end
    end
    return val
end


function show_positive_taxes(taxes::Array{Float64, 2})
    print("positive taxes: ")
    n, m = size(taxes)
    for i in 1:n
        for j in 1:m
            if taxes[i, j] > 0.01
                print("t[", i, ", ", j, "] = ", taxes[i, j], "  ")
            end
        end
    end
    print("\n")
end


function increase_taxes(instance::Instance, taxes::Array{Float64, 2})::Array{Float64, 2}
    for arc in shuffle!(collect(keys(instance.A1)))
        # For each arc in A1, given the other taxes,
        # we maximise the benefit made by this arc
        
        n1, n2 = arc
        initial_tax = taxes[n1, n2]
        current_tax = initial_tax
        
        previous_obj = 0 
        obj = compute_obj_value(instance, taxes, false)

        # we get all the (o, d) that go through this arc in a shortest path, 
        # and we get the max tax that this o-d allows before the clients change path 
        paths_through_a = compute_paths_through_a(instance, arc, taxes)
        possible_taxes = sort!(collect(values(paths_through_a)))
        iter = 1
        while obj >= previous_obj && iter <= length(possible_taxes)
            # while we improve the objective, we raise the tax
            previous_tax = current_tax
            current_tax = initial_tax + possible_taxes[iter] - EPS 
            # we remove epsilon to make sure that given 2 paths with same cost, 
            # cliens will choose paths with taxes
            if current_tax > 0 && abs(current_tax - previous_tax) > 10*EPS
                previous_obj = obj
                taxes[n1, n2] = current_tax
                obj = compute_obj_value(instance, taxes, false)
                if obj < previous_obj
                    # undo changes, because they worsened the objective
                    taxes[n1, n2] = previous_tax
                else
                    println("raising tax on arc ($n1, $n2) from $(round(Int, previous_tax)) to $(round(Int, current_tax))")
                end
            end
            iter += 1
        end
    end
    return taxes
end


function clean_taxes(taxes::Array{Float64, 2})
    n, m = size(taxes)
    for i in 1:n
        for j in 1:m
            taxes[i, j] = round(taxes[i, j])
        end
    end
end


function heuristic(instance::Instance, verbose::Bool=false)
    previous_val = 0
    n = instance.n
    nb_start = 5
    nb_iter = 5
    best_result = 0
    best_taxes = zeros(n, n)
    for i in 1:nb_start
        println("Start n°", i)
        taxes = zeros(n, n)
        current_val = 0
        for j in 1:nb_iter
            taxes = increase_taxes(instance, taxes)
            current_val = compute_obj_value(test_instance, taxes, (j==nb_iter && verbose))
        end
        if current_val > best_result
            best_result = current_val
            best_taxes = copy(taxes)
        end
    end

    # Undo the effect of epsilon
    clean_taxes(best_taxes)
    best_result = round(best_result)
    clean_taxes(best_taxes)
    best_result = round(best_result)

    return best_result, best_taxes
    return best_result, best_taxes
end