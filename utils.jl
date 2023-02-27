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
