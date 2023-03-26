# Utility functions

include("instance.jl")

function get_delta_plus_A1(instance::Instance, v::Int)
    delta_plus = Tuple{Int,Int}[]
    for (i,j) in eachindex(instance.A1)
        if i == v
            push!(delta_plus, (i,j))            
        end
    end
    return delta_plus
end

function get_delta_plus_A2(instance::Instance, v::Int)
    delta_plus = Tuple{Int,Int}[]
    for (i,j) in eachindex(instance.A2)
        if i == v
            push!(delta_plus, (i,j))            
        end
    end
    return delta_plus
end

function get_delta_minus_A1(instance::Instance, v::Int)
    delta_minus = Tuple{Int,Int}[]
    for (i,j) in eachindex(instance.A1)
        if j == v
            push!(delta_minus, (i,j))            
        end
    end
    return delta_minus
end

function get_delta_minus_A2(instance::Instance, v::Int)
    delta_minus = Tuple{Int,Int}[]
    for (i,j) in eachindex(instance.A2)
        if j == v
            push!(delta_minus, (i,j))            
        end
    end
    return delta_minus
end


function compute_obj_value(instance::Instance, taxes::Array{Float64, 2}, show_result::Bool=false)
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


function compute_paths_through_a(instance::Instance, arc::Tuple{Int, Int}, taxes::Array{Float64, 2})::Dict{}
    # returns a dictionnary with couples (origin, destination) as keys, and the max tax for (ori, des) on the arc as value
    n = instance.n
    adj_mat = build_dist_mat(instance.n, instance.A1, instance.A2, taxes)
    # ad_mat_without_taxed_arcs = build_dist_mat(instance.n, instance.A1, instance.A2, INF .* ones(n, n))
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
                
                cost, adj_mat[n1, n2] = adj_mat[n1, n2], INF
                next_dist = dijkstra_shortest_paths(instance.g, ori, adj_mat).dists[des]
                adj_mat[n1, n2] = cost
                # next_dist = dijkstra_shortest_paths(instance.g, ori, ad_mat_without_taxed_arcs).dists[des]

                # the maximum tax is the gap between the cost of the second best path and the first one
                paths_through_a[(ori, des)] = next_dist - dist
            end
        end
    end
    return paths_through_a
end


function local_search(instance::Instance, start_taxes::Matrix, naive::Bool=false)
    taxable_arcs = collect(keys(instance.A1))
    taxes = start_taxes
    current_obj = compute_obj_value(instance, taxes, false)
    println("Start Profit from Taxes : $(current_obj)")
    improvement = 1
    iteration = 1
    tested = 1
    max_T_max = maximum(instance.T_max)

    # While we improve our objective
    while improvement > 0
        # println("\nIteration $iteration \n")
        # We choose the first random (strictly) improving move
        for a in shuffle(taxable_arcs)
            iteration += 1
            # For each taxable edge, a neighboor is the previous or next possible tax on this edge
            current_tax = round(taxes[a[1], a[2]])
            possible_taxes = Float64[]
            
            if naive
                # We just do +1 and -1 on the arc
                current_tax + 1 <= max_T_max ? push!(possible_taxes, current_tax + 1) : 0
                current_tax - 1 >= 0 ? push!(possible_taxes, current_tax - 1) : 0
            else
                # We compute the first inscreasing neighbor
                paths_through_a = compute_paths_through_a(instance, a, taxes)
                filtered_taxes = filter(item -> item > 1e-5, unique(collect(values(paths_through_a))))
                length(filtered_taxes) > 0 ? push!(possible_taxes, current_tax + sort!(filtered_taxes)[1]) : 0
                # We compute the first decreasing neighbor
                if current_tax > 0
                    taxes[a[1], a[2]] = 0
                    paths_through_a = compute_paths_through_a(instance, a, taxes)
                    possible_decreasing = sort!(filter(item -> item > 0, unique(collect(values(paths_through_a)))), rev=true)
                    first_decreasing = findfirst(item -> item < current_tax, possible_decreasing)
                    if first_decreasing !== nothing 
                        push!(possible_taxes, possible_decreasing[first_decreasing])
                        # println("Possible decreasing : $possible_decreasing")
                    end
                end
            end
            
            # If no neighbor, we continue
            if length(possible_taxes) < 1
                continue
            end
            # println("Current tax on $a : $current_tax")
            # println("Neighbor taxes on $a : $(round.(Int, possible_taxes))")

            # We now test both neighbors (if they exist)
            for neighbor in possible_taxes
                taxes[a[1], a[2]] = neighbor- EPS
                improvement = compute_obj_value(instance, taxes, false) - current_obj
                # println("Testing neighbor : tax $(neighbor) on arc $a with improvement $(improvement)")
                if improvement < 0
                    # Undo changes
                    taxes[a[1], a[2]] = current_tax
                else
                    taxes[a[1], a[2]] = round(taxes[a[1], a[2]])
                    # Update objective
                    current_obj += round(improvement)
                    if round(improvement) > 0
                        println("Current tax on $a : $current_tax")
                        println("Neighbor taxes on $a : $(round.(Int, possible_taxes))")
                        println("Improving neighbor : tax $(neighbor) on arc $a with improvement $(improvement)")
                        println("Obj is now $current_obj")
                    end
                    # No need to check the other neighbor
                    break
                end
                tested += 1
            end
            if improvement > 0
                # No need to check the other edges
                break                
            end
        end
        println()
    end
    # println("End Profit from Taxes : $(current_obj)")
    # println("Iterations : $iteration")
    # println("Tested neighbors : $tested")
    return taxes, current_obj
end