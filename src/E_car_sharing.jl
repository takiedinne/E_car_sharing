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


include("simulation/data_structure.jl")
include("simulation/Vars.jl")
include("simulation/util.jl")

include("simulation/E_carsharing_sim.jl")

include("Mixed_integer_program/Mixed_integer_programme.jl")


export E_carsharing_sim,
       Solution,
       generate_random_solution,
       is_feasible_solution,
       car_type,
       set_online_mode, 
       initialize_scenario,
       initialize_scenarios,
       get_stored_solution
end