module E_car_sharing

using ResumableFunctions
using SimJulia

using DataFrames
using CSV

using Graphs, MetaGraphs
using EzXML

using Distances
using StatsBase

using LinearAlgebra
using JuMP, Gurobi
using Tables

using DelimitedFiles
using Serialization
using BenchmarkTools

using SparseArrayKit
using GraphMakie
using GLMakie

include("simulation/data_structure.jl")
include("simulation/Vars.jl")
include("simulation/util.jl")

include("simulation/visualizations.jl")

include("simulation/E_carsharing_sim.jl")

include("Mixed_integer_program/Mixed_integer_programme.jl")

include("expriments.jl")


export E_carsharing_sim,
    Solution,
    generate_random_solution,
    is_feasible_solution,
    car_type,
    set_online_mode,
    initialize_scenario,
    initialize_scenarios,
    get_stored_solution,
    get_solutions,
    serve_requests_after_opening_station,
    set_walking_time,
    set_cost_factor,
    project_path,
    clean_up_cars_number!,
    solve_using_mixed_integer_program
end