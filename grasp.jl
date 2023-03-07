# Random solution generation and local search heuristic

include("heuristic.jl")

function generate_random_taxes(instance::Instance)
    # We generate different solutions and take the best one
    n = instance.n
    min_T_max = Int(minimum(instance.T_max))
    max_T_max = Int(maximum(instance.T_max))
    taxes = zeros(max_T_max - min_T_max, n, n)
    # For each taxable edge, we choose a random value
    for e in shuffle!(collect(keys(instance.A1)))
        for M in 1:(max_T_max - min_T_max)
            taxes[M, e[1], e[2]] = rand(1:M)
        end
    end
    obj = [compute_obj_value(instance, taxes[i, :, :]) for i in 1:(max_T_max - min_T_max)]
    best_taxes = argmax(obj)
    println("Profit from Taxes : $(obj) with best at $(obj[best_taxes])")
    return taxes[best_taxes, :, :], obj[best_taxes]
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
        println("\nIteration $iteration \n")
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
    println("End Profit from Taxes : $(current_obj)")
    println("Iterations : $iteration")
    println("Tested neighbors : $tested")
    return taxes, current_obj
end

function grasp(instance::Instance)
    # Get code from ECMA
end

# Tests
include("data/taxe_grille_7x11.txt")
test_instance = Instance(K, n_k, n, A_1, A_2)
println(test_instance)
random_taxes, obj = generate_random_taxes(test_instance)
show_positive_taxes(random_taxes)

taxes, obj = local_search(test_instance, deepcopy(random_taxes))
taxes, obj = local_search(test_instance, deepcopy(random_taxes), true)