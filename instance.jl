# Instance object and related functions

using Graphs, SimpleWeightedGraphs

const EPS = 1e-5
const INF = 1e6
const INF_TAX = 1e4

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
    g :: SimpleDiGraph{Int64}

    function Instance(K, n_k, n, A_1, A_2)
        k = size(K)[1]
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


function build_graph(n, dict_A1, dict_A2)
    adj_matrix = zeros(n, n)
    for (e, cost) in dict_A1
        adj_matrix[e[1], e[2]] = 1
    end
    for (e, cost) in dict_A2
        adj_matrix[e[1], e[2]] = 1
    end
    g = DiGraph(adj_matrix)
    return g
end


function build_dist_mat(n::Int64, dict_A1, dict_A2, taxes::Matrix) :: Matrix
    dist_matrix = INF .* ones(n, n)
    for (e, cost) in dict_A1
        dist_matrix[e[1], e[2]] = cost + taxes[e[1], e[2]]
    end
    for (e, cost) in dict_A2
        dist_matrix[e[1], e[2]] = cost
    end
    return dist_matrix
end


function compute_t_max(K, n_k, n, dict_A1, dict_A2, g)
    k = size(K)[1]
    T_max = zeros(k)
    inf_taxes = INF_TAX .* ones(n, n)

    adj_mat = build_dist_mat(n, dict_A1, dict_A2, zeros(n, n))
    adj_mat_taxes = build_dist_mat(n, dict_A1, dict_A2, inf_taxes)
    
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

# Tests 
# include("data/taxe_grille_7x11.txt")
# test_instance = Instance(K, n_k, n, A_1, A_2)
# println(test_instance)