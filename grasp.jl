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

function local_search(instance::Instance, start_taxes::Matrix)
    taxable_arcs = collect(keys(instance.A1))
    taxes = start_taxes
    current_obj = compute_obj_value(instance, taxes, false)
    println("Start Profit from Taxes : $(current_obj)")
    improvement = 1
    iteration = 1

    # While we improve our objective
    while improvement > 0
        # We choose the first random (strictly) improving move
        for a in shuffle(taxable_arcs)
            # For each taxable edge, a neighboor is the previous or next possible tax on this edge
            paths_through_a = compute_paths_through_a(instance, a, taxes)
            possible_taxes = sort!(unique(collect(values(paths_through_a))))
            # We get the position of our current tax in the possible taxes
            current_tax = taxes[e[1], e[2]]
            current_idx = findfirst(item -> item >= current_tax, possible_taxes)
            # We now test both neighbors (if they exist)
            neighbors = [current_idx + i for i in [-1,1] if (current_idx + i >= 1 && current_idx + i <= length(possible_taxes))]
            for neighbor in neighbors
                taxes[e[1], e[2]] = possible_taxes[neighbor]
                improvement = compute_obj_value(instance, taxes, false) - current_obj
                if improvement <= 0
                    # Undo changes
                    taxes[e[1], e[2]] = current_tax
                else
                    # Update objective
                    current_obj += improvement
                    # No need to check the other neighbor
                    break
                end
                iteration += 1
            end
            if improvement > 0
                # No need to check the other edges
                break                
            end
        end
    end
    println("End Profit from Taxes : $(current_obj)")
    println("Tested neighbors : $iteration")
    return taxes, current_obj
end

function grasp(instance::Instance)
    # Get code from ECMA
end

# Tests
include("data/taxe_grille_7x11.txt")
test_instance = Instance(K, n_k, n, A_1, A_2)
# for i in 1:10
#     random_taxes, obj = generate_random_taxes(test_instance)
# end
# show_positive_taxes(random_taxes)
