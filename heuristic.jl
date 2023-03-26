# Heuristic resolution

using JuMP
using CPLEX
using Random
using Statistics

include("instance.jl")
include("utils.jl")

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


function increase_taxes(instance::Instance, taxes::Array{Float64, 2}, verbose::Bool)::Array{Float64, 2}
    p = 1 # probability of keeping a tax that does not change the objective

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
                elseif obj - previous_obj <= EPS && rand() < 1 - p
                    # if the objective is the same, we keep changes with probability p
                    taxes[n1, n2] = previous_tax
                elseif verbose
                    println("raising tax on arc ($n1, $n2) from $(round(Int, previous_tax)) to $(round(Int, current_tax))")
                end
            end
            iter += 1
        end
    end
    return taxes
end

function generate_taxes(instance::Instance, nb_incr::Int, verbose::Bool=false)
    n = instance.n
    taxes = zeros(n, n)
    current_val = 0
    for j in 1:nb_incr
        taxes = increase_taxes(instance, taxes, verbose)
        current_val = compute_obj_value(instance, taxes, (j==nb_incr && verbose))
    end
    return taxes, current_val
end


function clean_taxes(taxes::Array{Float64, 2})
    n, m = size(taxes)
    for i in 1:n
        for j in 1:m
            taxes[i, j] = round(taxes[i, j])
        end
    end
end

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

function heuristic(instance::Instance, verbose::Bool=false)
    n = instance.n
    nb_start = 400
    nb_iter = 5
    best_result = 0
    best_taxes = zeros(n, n)
    all_taxes = zeros(nb_start)
    for i in 1:nb_start
        current_val = 0
        if rand() < 0.3
            taxes, current_val = generate_taxes(instance, nb_iter, verbose)
        else
            taxes, current_val = generate_random_taxes(instance)
        end
        taxes, current_val = local_search(instance, taxes)
        all_taxes[i] = current_val
        if current_val > best_result
            best_result = current_val
            best_taxes = copy(taxes)
        end
    end

    println("Moyenne des profits : ", mean(all_taxes), " et écart-type : ", std(all_taxes))

    # Undo the effect of epsilon
    clean_taxes(best_taxes)
    best_result = round(best_result)

    return best_result, best_taxes
end