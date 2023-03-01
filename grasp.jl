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
    tested = 1

    # While we improve our objective
    while improvement > 0
        println("\nIteration $iteration \n")
        # We choose the first random (strictly) improving move
        for a in shuffle(taxable_arcs)
            iteration += 1
            # For each taxable edge, a neighboor is the previous or next possible tax on this edge
            paths_through_a = compute_paths_through_a(instance, a, taxes)
            possible_taxes = sort!(unique(collect(values(paths_through_a))))
            if length(possible_taxes) < 1
                continue
            end
            println("Possible taxes on $a : $(round.(Int, possible_taxes))")
            # We get the position of our current tax in the possible taxes
            current_tax = round(taxes[a[1], a[2]])
            println("Current tax on $a : $current_tax")
            current_idx = findfirst(item -> item >= current_tax, possible_taxes)
            
            # We construct the neighbors
            neighbors = []
            if current_idx === nothing 
                # The current tax is superior to all values in possible taxes, the only neighbor is then the last value
                push!(neighbors, length(possible_taxes))
            else
                if cuurent_idx == 1
                    # The current tax is inferior to all values in possible taxes, the only neighbor is then the first value
                    push!(neighbors, 1)
                else
                    # There are possible values inferior and superior to the current tax
                    push!(neighbors, current_idx)
                    push!(neighbors, current_idx - 1)
                end
            end

            # We now test both neighbors (if they exist)
            # TODO : ADAPT FROM HERE
            for neighbor in neighbors
                taxes[a[1], a[2]] = possible_taxes[neighbor] - EPS
                improvement = compute_obj_value(instance, taxes, false) - current_obj
                println("Testing neighbor : tax $(possible_taxes[neighbor]) on arc $a with improvement $improvement")
                if improvement <= 0
                    # Undo changes
                    taxes[a[1], a[2]] = current_tax
                else
                    # Update objective
                    current_obj += improvement
                    println("Obj is now $current_obj")
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
# for i in 1:10
#     random_taxes, obj = generate_random_taxes(test_instance)
# end
# show_positive_taxes(random_taxes)

n = test_instance.n
taxes, obj = local_search(test_instance, zeros(n, n))