time_slot_length = 5 # min 
scenario_path_list = ["Data/Scenarios_1000_greaterthan2/Output1000_1.txt"#= ,
    "Data/Scenarios_1000_greaterthan2/Output1000_2.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_3.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_4.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_5.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_6.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_7.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_8.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_9.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_10.txt" =#]

J = get_potential_locations() #set of potentiel station
stations_idx = Dict((J[i] => i) for i in eachindex(J)) # useful to get the index of station inside J vector
f_j = [get_prop(manhaten_city_graph, i, :charging_station_base_cost) for i in J] # fixed cost of each station
g = vehicle_specific_values[Smart_ED][:car_cost] # the cost of purchasing the car 
C_j = [get_prop(manhaten_city_graph, j, :max_number_of_charging_points) for j in J] # the capacity of each station
T = collect(Int64, 1:300) # set of time slots
S = collect(1:1) # we are  only considering one scenario at a time
q_s = [1] # the probabilty of each scenario here we are only consedring  one scenario

β = vehicle_specific_values[Smart_ED][:battery_capacity] # battery capacity
δ_ij = hcat([[get_trip_battery_consumption(J[i], J[j], Smart_ED) * β / 100 for i in eachindex(J)] for j in eachindex(J)]...) # batteru usage between the stations

d_ij = hcat([[get_trip_duration(J[i], J[j]) for i in eachindex(J)] for j in eachindex(J)]...) # driving time between stations
#d_w = hcat([[get_walking_time(i, J[j]) for i in 1:nv(manhaten_city_graph)] for j in eachindex(J)]...)# walking time between each node and each station
#CSV.write("walking_time.csv", Tables.table(d_w), writeheader=false)
dw_ij = readdlm(project_path("Data/MIP/walking_time.csv"), ',', Float64)
β_w = 5 # maximum allowed walking time

charging_rate = vehicle_specific_values[Smart_ED][:fast_charging_rate]

solution_list = Array{Solution,1}()
objective_value_list = Array{Float64,1}()

function create_solution_for_simulation(sol_df::DataFrame)
    sol_y = filter(row -> occursin("y", row.x), sol_df)
    sol_u_h = filter(row -> occursin("u_h", row.x), sol_df)
    sol_L_0 = filter(row -> occursin("L_0", row.x), sol_df)
    
    return Solution( convert.(Bool, sol_y.val),  convert.(Int, sol_L_0.val), [convert.(Bool, sol_u_h.val)])
end

function get_solutions()
    scenario_counter = 1
    mip = Model()
    for path in scenario_path_list
        println("sceario [$scenario_counter]: start working with scenario ")
        scenario_list = [requests_as_dataframe(path)]
        for i in eachindex(scenario_list)
            scenario_list[i].scenarioId = i .* ones(Int, nrow(scenario_list[i]))
        end

        K_s = scenario_list
        K = vcat(scenario_list...) #set of all request in all scenarios
        Δ_k = [get_trip_duration(req.ON, req.DN) * 1.1 for req in eachrow(K)] #time threshold of the trip duration for each request

        Hs = [get_feasible_paths(scenario, J, β_w) for scenario in scenario_list]# all feasible paths
        #add scenario id 
        for i in eachindex(Hs)
            Hs[i].scenario_id = i .* ones(Int, nrow(Hs[i]))
        end
        H = vcat(Hs...)


        #generate the parameters
        b = zeros(Bool, nrow(H), length(J), length(T))
        μ = zeros(Bool, nrow(H), length(J), length(T))
        λ = zeros(Bool, nrow(H), length(J), length(T))

        for h in 1:nrow(H)
            # set bata
            start_station_id = stations_idx[H.origin_station[h]]
            destination_station_id = stations_idx[H.destination_station[h]]
            starting_time_slot = H.start_driving_time[h]
            b[h, start_station_id, starting_time_slot] = 1

            # set μ
            start_charging_time_slot = H.arriving_time[h]
            end_charging_time_slot = H.arriving_time[h] + ceil(Int, (δ_ij[start_station_id, destination_station_id] / charging_rate) / 60 / time_slot_length)
            μ[h, destination_station_id, start_charging_time_slot:end_charging_time_slot] .= 1

            # set λ
            λ[h, destination_station_id, end_charging_time_slot+1] = 1
        end

        # begin mixed integer programme
        mip = Model()

        println("sceario [$scenario_counter]: constructing the mixed integer program ...")
        # define the variables
        if isfile("Data/MIP/programs_file/E_carsharing_mip_scenario_$scenario_counter.mof.json")
            mip = read_from_file("Data/MIP/programs_file/E_carsharing_mip_scenario_$scenario_counter.mof.json")
        else
            @variable(mip, u_h[1:nrow(H)], Bin) # u_h 1 if the trip in H is used 0 otherwise
            @variable(mip, L_ts[eachindex(J), eachindex(T), eachindex(S)] >= 0, Int)
            @variable(mip, L_0[eachindex(J)], Int)
            @variable(mip, y[eachindex(J)], Bin)
        
            # List constraints
        
            #constraints (2)
            #@constraint(mip, c2[i=1:nrow(K)], sum(u_h[j] for j in findall(x -> x == K.reqId[i], H.req)) <= 1 )
            @constraint(mip, c2[i=1:nrow(K)], sum(u_h[j] for j in findall((H.req .== K.reqId[i]) .& (H.scenario_id .== K.scenarioId[i]))) <= 1)
        
            #constraints (3)
            @constraint(mip, c3_1[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.origin_station[i]]])
            @constraint(mip, c3_2[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.destination_station[i]]])
        
            #constraint (4)
            @constraint(mip, c4[j=eachindex(J), s=eachindex(S), t=eachindex(T)], sum(b[h, j, t] * u_h[h] for h in findall(H.scenario_id .== s)) <= L_ts[j, t, s])
        
            #constraint (5)
            @constraint(mip, c5[j=eachindex(J), s=eachindex(S), t=eachindex(T)],
                L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h]
                                    for h in findall(H.scenario_id .== s)) <= C_j[j] * y[j])
        
            #constraint (6)
            @constraint(mip, c6[j=eachindex(J), s=eachindex(S), t=2:length(T)], L_ts[j, t, s] == L_ts[j, t-1, s] + sum((λ[h, j, t] - b[h, j, t-1]) * u_h[h] for h in findall(H.scenario_id .== s)))
        
            #constraint (7)
            @constraint(mip, c7[j=eachindex(J), s=eachindex(scenario_list)], L_ts[j, 1, s] == L_0[j])
        
            #constraint (8)
            @constraint(mip, c8_1[j=eachindex(J), t=eachindex(T), s=eachindex(scenario_list)], L_ts[j, t, s] <= C_j[j] * y[j])
            @constraint(mip, c8_2[j=eachindex(J), t=eachindex(T), s=eachindex(scenario_list)], L_ts[j, t, s] >= 0)
        
            #constraint (9)
            @constraint(mip, c9[j=eachindex(J)], L_0[j] >= 0)
        
            #constraint (10) and (11) are precised whene defining the variables
        
        
            # the objective function
            @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in findall(H.scenario_id .== s)) for s in eachindex(S)) -
                                 sum(f_j[j] * y[j] for j in eachindex(J)) / cost_factor -
                                 g * sum(L_0[j] for j in eachindex(J)) / cost_factor)
        
        
            #save the programme
            write_to_file(mip, "Data/MIP/programs_file/E_carsharing_mip_scenario_$scenario_counter.mof.json")
        end
        
        set_optimizer(mip, Gurobi.Optimizer)
        set_time_limit_sec(mip, 1200.0)

        println("sceario [$scenario_counter]: solving the mixed integer program")
        optimize!(mip)

        obj_val = objective_value(mip)

        sol_df = DataFrame(x=name.(all_variables(mip)), val=JuMP.value.(all_variables(mip)))

        sol = create_solution_for_simulation(sol_df)
        scenario_counter += 1
        
        push!(solution_list, sol)
        push!(objective_value_list, objective_value(mip))
    end

    println("Saving the solutions")
    for i in eachindex(solution_list)
        serialize("Data/MIP/solutions/solution_$i.jls", solution_list[i])
    end

    serialize("Data/MIP/solutions/Objective_values.jls", objective_value_list)

    println("End traitement pass to the validation")

    mip
end


