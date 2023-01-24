# Instance object and related functions

include("data/taxe_grille_2x3.txt")

struct Instance
    # Commodity parameters
    K :: Int
    origins :: Vector{Int}
    destinations :: Vector{Int}
    demands :: Vector{Int}

    # Graph parameters
    n :: Int
    n_A1 :: Int
    A1 :: Dict{Tuple{Int, Int}, Float64}
    n_A2 :: Int
    A2 :: Dict{Tuple{Int, Int}, Float64}
    T_max :: Float64

    function Instance(K, n_k, n, A_1, A_2)
        k = size(K)[2]
        origins = K[:,1]
        destinations = K[:,2]
        dict_A1 = matrix_to_dict(A_1)
        dict_A2 = matrix_to_dict(A_2)
        T_max = maximum(values(dict_A2)) - minimum(values(dict_A1))
        new(k, origins, destinations, n_k, n, length(dict_A1), dict_A1, length(dict_A2), dict_A2, T_max)
    end
end

function Base.show(io::IO, instance::Instance)
    println("Instance : ")
    println(string("K = ", instance.K))
    ODD = zeros(Int, instance.K, 3)
    ODD[:,1] = instance.origins
    ODD[:,2] = instance.destinations
    ODD[:,3] = instance.demands
    println("Origin-Destination-Demand matrix : ")
    display(ODD)
    println()
    println(string("n = ", instance.n))
    println(string("n_A1 = ", instance.n_A1))
    println("A1 : ")
    display(instance.A1)
    println()
    println("A2 : ")
    display(instance.A2)
    println()
    println(string("T_max = ", instance.T_max))
    println()
end

function matrix_to_dict(matrix::Matrix{Int})
    dict = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:length(matrix[:,1])
        i, j, k = matrix[i,:]
        dict[i,j] = k
    end
    return dict
end

# Test
test_instance = Instance(K, n_k, n, A_1, A_2)
println(test_instance)
