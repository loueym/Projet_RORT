# Heuristic resolution

using JuMP
using CPLEX

include("instance.jl")

function compute_shortest_path(instance::Instance, journey_index::Int64, taxes)
    # Model
    model = Model(CPLEX.Optimizer) 
    
    # Variables
    @variable(model, x1[(i,j) in eachindex(instance.A1)] >= 0)
    @variable(model, x2[(i,j) in eachindex(instance.A2)] >= 0)

    # add constraints

    # Objective
    @objective(
        model, 
        Min, 
        (
            sum(x1[i, j]*(instance.A1[i, j] + get(taxes, (i, j), 0)) for (i, j) in keys(instance.A1)) 
            + sum(x2[i, j]*instance.A2[i, j] for (i, j) in keys(instance.A2))
        )
    )

end
