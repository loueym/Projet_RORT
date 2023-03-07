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