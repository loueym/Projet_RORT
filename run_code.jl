include("data/taxe_grille_7x11.txt")
include("heuristic.jl")

test_instance = Instance(K, n_k, n, A_1, A_2)
# println(test_instance)

# g = test_instance.g
# for i in 1:n
#     for j in i+1:n
#         println(i, " -> ", j, " ", get_weight(g, i, j))
#         println(j, " -> ", i, " ", get_weight(g, i, j))
#     end
# end

# dijk = dijkstra_shortest_paths(g, [1, 2], build_adj_mat(n, test_instance.A1, test_instance.A2))
# println(dijk.dists)
# dijk_enum = enumerate_paths(dijkstra_shortest_paths(g, 2, build_adj_mat(n, test_instance.A1, test_instance.A2)), 8)
# println(dijk_enum)

value, taxes = heuristic(test_instance)
println("Value after heuristic: ", value)
show_positive_taxes(taxes)