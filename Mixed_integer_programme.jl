using LinearAlgebra
using JuMP
using Gurobi
using GLPK
using DataFrames
using CSV
using Tables
using DelimitedFiles
using BenchmarkTools
using Serialization

include("E_carsharing_sim.jl")

function create_solution_for_simulation(sol_df::DataFrame)
    sol_y = filter(row->occursin("y", row.x), sol_df)
    sol_u_h = filter(row->occursin("u_h", row.x), sol_df)
    sol_L_0 = filter(row->occursin("L_0", row.x), sol_df)

    # solution fields
    open_stations_ids = J[findall(val -> val == 1, sol_y.val)]
    initial_car_number = sol_L_0.val[findall(val -> val == 1, sol_y.val)]
    selected_paths = sol_u_h.val

    return Solution(open_stations_ids, initial_car_number, selected_paths)
end

time_slot_length = 5 # min 
scenario_path_list = [  "Data/Scenarios_1000_greaterthan2/Output1000_1.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_2.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_3.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_4.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_5.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_6.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_7.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_8.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_9.txt",
                        "Data/Scenarios_1000_greaterthan2/Output1000_10.txt"]



J = get_potential_locations() #set of potentiel station
f_j = [get_prop(manhaten_city_graph, i, :charging_station_base_cost) for i in J] # fixed cost of each station
g = vehicle_specific_values[Smart_ED][:car_cost] # the cost of purchasing the car 
C_j = [get_prop(manhaten_city_graph, j, :max_number_of_charging_points) for j in J] # the capacity of each station
T = collect(Integer, 1:24*60/5+10) # set of time slots
S = collect(1:1)
q_s = [1] # the probabilty of each scenario here we are consedring only one scenario
δ_ij = hcat([[get_trip_battery_consumption(J[i], J[j], Smart_ED) for i in 1:length(J)] for j in 1:length(J)]...) # batteru usage between the stations
β = vehicle_specific_values[Smart_ED][:battery_capacity] # battery capacity
d_ij = hcat([[get_trip_duration(J[i], J[j]) for i in 1:length(J)] for j in 1:length(J)]...) # driving time between stations
#d_w = hcat([[get_walking_time(i, J[j]) for i in 1:nv(manhaten_city_graph)] for j in 1:length(J)]...)# walking time between each node and each station
#CSV.write("walking_time.csv", Tables.table(d_w), writeheader=false)
dw_ij = readdlm("Data/MIP/walking_time.csv", ',', Float64)
β_w = 5 # maximum allowed walking time
  
charging_rate = vehicle_specific_values[Smart_ED][:fast_charging_rate]

scenario_counter = 1

solution_list = Array{Solution, 1}()
objective_value_list = Array{Float64, 1}()

for path in scenario_path_list
    
    println("sceario [$i]: start working with scenario ")
    scenario_list = [scenario_as_dataframe(path)]
    for i in 1:length(scenario_list)
        scenario_list[i].scenarioId =  i .* ones(Int, nrow(scenario_list[i]))
    end

    K_s = scenario_list
    K = vcat(scenario_list...) #set of all request in all scenarios
    Δ_k = [get_trip_duration(req.ON, req.DN) * 1.1 for req in eachrow(K)] #time threshold of the trip duration for each request

    Hs = [get_all_requests_feasible_paths(scenario, J, β_w) for scenario in scenario_list]# all feasible paths
    #add scenario id 
    for i in 1:length(Hs)
        Hs[i].scenario_id = i .* ones(Int, nrow(Hs[i]))
    end
    H =  vcat(Hs...)

    #generate the parameters
    b = zeros(Bool, nrow(H), length(J), length(T))
    μ = zeros(Bool, nrow(H), length(J), length(T))
    λ = zeros(Bool, nrow(H), length(J), length(T))
    
    for h in 1:nrow(H)
        # set bata
        start_station_id = findfirst(x -> x == H.origin_station[h], J)
        destination_station_id = findfirst(x -> x == H.destination_station[h], J)
        starting_time_slot = floor(Int, H.start_driving_time[h] / 5) + 1
        b[h, start_station_id, starting_time_slot] = 1
        
        # set μ
        start_charging_time_slot = floor(Int, H.arriving_time[h] / 5) + 1
        end_charging_time_slot =  floor(Int, (H.arriving_time[h] +  δ_ij[start_station_id, destination_station_id] / charging_rate ) / 5) + 1 
        μ[h, destination_station_id, start_charging_time_slot:end_charging_time_slot] .= 1

        # set λ
        λ[h, destination_station_id, end_charging_time_slot + 1 ] = 1
    end 
    
    # begin mixed integer programme
    mip = Model()
    
    println("sceario [$i]: constructing the mixed integer program ...")
    # define the variables 
    @variable(mip, u_h[1:nrow(H)], Bin) # u_h 1 if the trip in H is used 0 otherwise
    @variable(mip, L_ts[1:length(J), 1:length(T), 1:length(S)] >= 0, Int) 
    @variable(mip, L_0[1:length(J)], Int)
    @variable(mip, y[1:length(J)], Bin)
    
    # List constraints

    #constraints (2)
    @constraint(mip, c2[i=1:nrow(K)], sum(u_h[j] for j in findall(x -> x == K.reqId[i], H.req)) <= 1 )

    #constraints (3)
    @constraint(mip, c3_1[i=1:nrow(H)], u_h[i] <= y[findfirst(x-> x == H.origin_station[i], J)] )
    @constraint(mip, c3_2[i=1:nrow(H)], u_h[i] <= y[findfirst(x-> x == H.destination_station[i], J)] )

    #constraint (4)
    @constraint(mip, c4[j=1:length(J), s=1:length(S), t=1:length(T)], 
                            sum(b[h, j, t] * u_h[h] for h in findall(x-> x == s, H.scenario_id)) <= L_ts[j, t, s] )

    #constraint (5)
    @constraint(mip, c5[j = 1:length(J), s=1:length(S), t=1:length(T)],
                        L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h] 
                        for h in findall(x-> x == s, H.scenario_id)) <= C_j[j] * y[j])

    #constraint (6)
    @constraint(mip, c6[j = 1:length(J), s=1:length(S), t=2:length(T)], L_ts[j, t, s] == L_ts[j, t-1, s] +  sum((λ[h, j, t] - b[h, j, t-1]) * u_h[h] for h in findall(x-> x == s, H.scenario_id)) )

    #constraint (7)
    @constraint(mip, c7[j=1:length(J), s=1:length(scenario_list)], L_ts[j, 1, s] == L_0[j])

    #constraint (8)
    @constraint(mip, c8_1[j=1:length(J), t=1:length(T), s=1:length(scenario_list)], L_ts[j, t, s] <= C_j[j] * y[j])
    @constraint(mip, c8_2[j=1:length(J), t=1:length(T), s=1:length(scenario_list)], L_ts[j, t, s] >= 0)

    #constraint (9)
    @constraint(mip, c9[j=1:length(J)], L_0[j] >= 0)

    #constraint (10) and (11) are precised whene defining the variables


    # the objective function
    @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in 1:nrow(H))  for s in 1:length(S)) - 
                sum( f_j[j] * y[j] for j in 1:length(J)) / cost_factor - 
                g * sum(L_0[j] for j in 1:length(J)) / cost_factor )


    #save the programme
    write_to_file(mip, "Data/MIP/E_carsharing_mip_scenario_$scenario_counter.mof.json")
    
    #mip = read_from_file("Data/MIP/E_carsharing_mip_scenario1.mof.json")
    set_optimizer(mip, Gurobi.Optimizer)
    set_time_limit_sec(mip, 1200.0)

    println("sceario [$i]: solving the mixed integer program")
    optimize!(mip)

    obj_val = objective_value(mip)

    sol_df = DataFrame(x = name.(all_variables(mip)), val = JuMP.value.(all_variables(mip)) )

    sol = create_solution_for_simulation(sol_df)
    global scenario_counter += 1

    push!(solution_list, sol)
    push!(objective_value_list, objective_value(mip))
end

println("Saving the solutions")
for i in 1:length(solution_list)
    serialize("Data/MIP/solutions/solution_$i.jls", solution_list[i])
end

serialize("Data/MIP/solutions/Objective_values.jls", objective_value_list)

println("End traitement pass to the validation")