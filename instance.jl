# Instance object and related functions

using LightGraphs, SimpleWeightedGraphs

struct Instance
    # Commodity parameters
    K :: Int
    origins :: Vector{Int}
    destinations :: Vector{Int}
    demands :: Vector{Int}
    T_max :: Vector{Float64}

    # Graph parameters
    n :: Int
    n_A1 :: Int
    A1 :: Dict{Tuple{Int, Int}, Float64}
    n_A2 :: Int
    A2 :: Dict{Tuple{Int, Int}, Float64}
    g :: SimpleWeightedDiGraph{Int64,Float64}

    function Instance(K, n_k, n, A_1, A_2)
        k = size(K)[2]
        origins = K[:,1]
        destinations = K[:,2]
        dict_A1 = matrix_to_dict(A_1)
        dict_A2 = matrix_to_dict(A_2)
        g = build_graph(n, dict_A1, dict_A2)
        T_max = compute_t_max(K, n_k, n, dict_A1, dict_A2, g)
        new(k, origins, destinations, n_k, T_max, n, length(dict_A1), dict_A1, length(dict_A2), dict_A2, g)
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

function inf_mat(mat::Array{Float64, 2}, inf::Union{Float64,Int64})
    n, m = size(mat)
    for i in 1:n
        for j in 1:m
            mat[i, j] = inf
        end
    end
end


function build_graph(n, dict_A1, dict_A2)
    adj_matrix = zeros(n, n)
    eps = 0.0001
    inf = 1000000
    inf_mat(adj_matrix, inf)
    for (e, cost) in dict_A1
        if cost == 0
            cost = eps
        end
        adj_matrix[e[1], e[2]] = cost
    end
    for (e, cost) in dict_A2
        if cost == 0
            cost = eps
        end
        adj_matrix[e[1], e[2]] = cost
    end
    g = SimpleWeightedDiGraph(adj_matrix)
    for i in 1:n
        for j in i+1:n
            if adj_matrix[i, j] == inf
                rem_edge!(g, i, j)
            end
            if adj_matrix[j, i] == 1000000
                rem_edge!(g, j, i)
            end
        end
    end
    return g
end


function build_adj_mat(n::Int64, dict_A1, dict_A2, taxes::Array{Float64, 2})::Array{Float64, 2}
    adj_matrix = zeros(n, n)
    inf = 1000000
    inf_mat(adj_matrix, inf)
    for (e, cost) in dict_A1
        adj_matrix[e[1], e[2]] = cost + taxes[e[1], e[2]]
    end
    for (e, cost) in dict_A2
        adj_matrix[e[1], e[2]] = cost
    end
    return adj_matrix
end


function compute_t_max(K, n_k, n, dict_A1, dict_A2, g)
    T_max = Vector{Float64}([0 for i in 1:n])
    inf_taxes = zeros(n, n)
    inf_mat(inf_taxes, 10000)
    adj_mat = build_adj_mat(n, dict_A1, dict_A2, zeros(n, n))
    adj_mat_taxes = build_adj_mat(n, dict_A1, dict_A2, inf_taxes)
    k = size(K)[2]
    for i in 1:k
        ori = K[i,1]
        des = K[i,2]
        demand = n_k[i]
        dist = dijkstra_shortest_paths(g, ori, adj_mat).dists[des]
        dist_taxes = dijkstra_shortest_paths(g, ori, adj_mat_taxes).dists[des]
        T_max[i] = dist_taxes - dist
    end
    return T_max
end