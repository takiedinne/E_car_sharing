J = get_potential_locations() #set of potentiel station
stations_idx = Dict((J[i] => i) for i in eachindex(J)) # useful to get the index of station inside J vector
f_j = [get_prop(manhaten_city_graph, i, :charging_station_base_cost) for i in J] # fixed cost of each station
g = vehicle_specific_values[Smart_ED][:car_cost] # the cost of purchasing the car 
C_j = [get_prop(manhaten_city_graph, j, :max_number_of_charging_points) for j in J] # the capacity of each station
T = collect(Int64, 1:300) # set of time slots
S = collect(1:1) # to redefine
q_s = [1] # to redefine

β = vehicle_specific_values[Smart_ED][:battery_capacity] # battery capacity
# battery usage between the stations. get_trip_duration return percentage so we need to converted again to watt
δ_ij = hcat([[get_trip_battery_consumption(J[i], J[j], Smart_ED) * β / 100 for i in eachindex(J)] for j in eachindex(J)]...) 

#d_ij = hcat([[get_trip_duration(J[i], J[j]) / time_slot_length for i in eachindex(J)] for j in eachindex(J)]...) # driving time between stations
#dw_ij = hcat([[get_walking_time(i, J[j]) for i in 1:nv(manhaten_city_graph)] for j in eachindex(J)]...)# walking time between each node and each station
#CSV.write("walking_time.csv", Tables.table(d_w), writeheader=false)
#dw_ij = readdlm(project_path("Data/MIP/walking_time.csv"), ',', Float64)
β_w = 5 # maximum allowed walking time

H = DataFrame() #set all feasible paths
sc_fps_range = [] # useful to filter on the h of each scenario
charging_rate = vehicle_specific_values[Smart_ED][:fast_charging_rate]

solution_list = Array{Solution,1}()
objective_value_list = Array{Float64,1}()
"""
    create_solution_for_simulation(sol_df::DataFrame, H::DataFrame)

create a solution object from the solution of the MIP after optimization, where the name of variablea
and their values are stored in `H`` dataframe

"""
function create_solution_for_simulation(sol_df::DataFrame, H::DataFrame)
    sol_y = filter(row -> occursin("y", row.x), sol_df)
    sol_u_h = filter(row -> occursin("u_h", row.x), sol_df)
    sol_L_0 = filter(row -> occursin("L_0", row.x), sol_df)
    
    #selected paths  special traitement
    scenarios = unique(H.scenarioId) 

    selected_paths = [convert.(Bool, round.(sol_u_h.val[findall(x->x == i, H.scenarioId)])  ) for i in scenarios]
    
    return Solution(convert.(Bool, round.(sol_y.val)),  convert.(Int, round.(sol_L_0.val)), selected_paths)
end

""" 
    solve the MIP for the given scenarios
    @scenarios: is the list for the scenarios after they have been initialized
    @genarated_scenarios (optional): if the scenarios are amoung the generated scenario or false when we are using the already exist scenarios
"""
function solve_using_mixed_integer_program(scenarios::Vector{Scenario}; genarated_scenarios = false, mip_file_path="")
    
    # we make a convention that the file is names as "E_carsharing_mip_1_to_$(nbr_of_scenarios).mof.json"
    # so is we wanted to solve 5 scenarions so the name will be  "E_carsharing_mip_1_to_5.mof.json"
    if mip_file_path == ""
        mip_file_suffix = genarated_scenarios ? "_generated_$(length(scenarios))" : "1_to_$(length(scenarios))"
        mip_file_path = "Data/MIP/programs_file/E_carsharing_mip_$(mip_file_suffix)_$(maximum_walking_time)_walking_time.mof.json"
    end
    
    # The mixed integer program 
    mip = Model()

    # we need H variable in both case so we define it here
    # add column to the different scenarios to store the scenario id
    for i in eachindex(scenarios)
        scenarios[i].request_list.scenarioId = i .* ones(Int, nrow(scenarios[i].request_list)) 
        scenarios[i].feasible_paths.scenarioId = i .* ones(Int, nrow(scenarios[i].feasible_paths)) 
    end
    
    #Hs = [scenarios[i].feasible_paths  for i in eachindex(scenarios)]# all feasible paths
    #H = vcat(Hs...)
    H = vcat([scenarios[i].feasible_paths  for i in eachindex(scenarios)]...)
    
    # feasible paths range (to avoid filtring at each time)
    sc_fps_range= Array{UnitRange, 1}(undef, length(scenarios))
    curent_index = 1
    for sc in eachindex(scenarios)
        #get the range of feasible paths for the corresponding scenario
        sc_fps_range[sc] = UnitRange(curent_index, curent_index + nrow(scenarios[sc].feasible_paths) - 1)
        curent_index += nrow(scenarios[sc].feasible_paths)
    end

    if isfile(mip_file_path)
        mip = read_from_file(mip_file_path)
    else
        
        @info "[Solving the MIP]: creating the MIP ..."
        #set the different variables for the MIP and try to preserve the same names as in the paper
        #K_s = [scenarios[i].request_list for i in eachindex(scenarios)]
        #K = vcat(K_s...) #set of all request in all scenarios
        K = vcat([scenarios[i].request_list for i in eachindex(scenarios)]...)
        #Δ_k = [get_trip_duration(req.ON, req.DN) / time_slot_length * 1.1 for req in eachrow(K)] #time threshold of the trip duration for each request
        
        S = collect(1:length(scenarios)) 
        # at this stage we are not taking in consideration the probability of each scenario
        q_s = 1 / length(scenarios) .* ones(length(scenarios)) 

        @info "[Solving the MIP]: creating the variables ..."
        #generate the parameters
        b = SparseArray{Bool}(undef, nrow(H), length(J), length(T)) 
        μ = SparseArray{Bool}(undef, nrow(H), length(J), length(T)) 
        λ = SparseArray{Bool}(undef, nrow(H), length(J), length(T))
        
        for h in 1:nrow(H)
            # set bata
            start_station_id = stations_idx[H.origin_station[h]]
            destination_station_id = stations_idx[H.destination_station[h]]
            starting_time_slot = H.start_driving_time[h]
            b[h, start_station_id, starting_time_slot] = 1
            #@benchmark b[h, stations_idx[H.origin_station[h]], H.start_driving_time[h]] = 1
            
            # set μ
            start_charging_time_slot = H.arriving_time[h]
            end_charging_time_slot = H.arriving_time[h] + ceil(Int, 
                                (δ_ij[start_station_id, destination_station_id] / charging_rate) / 60 / time_slot_length)
            μ[h, destination_station_id, start_charging_time_slot:end_charging_time_slot] .= 1
            
            # set λ
            λ[h, destination_station_id, end_charging_time_slot+1] = 1
        end
        
        # define the variables
        @variable(mip, u_h[1:nrow(H)], Bin) # u_h 1 if the trip in H is used 0 otherwise
        @variable(mip, L_ts[eachindex(J), eachindex(T), eachindex(S)] >= 0, Int)
        @variable(mip, L_0[eachindex(J)] >=0 , Int)
        @variable(mip, y[eachindex(J)], Bin)
        
        # List constraints
        @info "[Solving the MIP]: creating the constraints ..."
        #constraints (2)
        #= @constraint(mip, c2[i=1:nrow(K)], 
                    sum(u_h[j] for j in findall((H.req .== K.reqId[i]) 
                    .& 
                    (H.scenarioId .== K.scenarioId[i]))) <= 1) =#

        @constraint(mip, c2[i=1:nrow(K)], sum(u_h[j] for j in (K.fp[i] .+ (1000 * (K.scenarioId[i] - 1 ))) ) <= 1)

        #constraints (3)
        @constraint(mip, c33_1[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.origin_station[i]]])
        @constraint(mip, c33_2[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.destination_station[i]]])
        
        #constraint (4)
       #=  @constraint(mip, c4[j=eachindex(J), s=eachindex(S), t=eachindex(T)], 
                        sum(b[h, j, t] * u_h[h] for h in findall(H.scenarioId .== s)) <= L_ts[j, t, s]) =#
                 
        @constraint(mip, c4[j=eachindex(J), s=eachindex(S), t=eachindex(T)], 
                        sum(b[h, j, t] * u_h[h] for h in sc_fps_range[s]) <= L_ts[j, t, s])
              
        #constraint (5)
        #= @constraint(mip, c5[j=eachindex(J), s=eachindex(S), t=eachindex(T)],
            L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h]
                                for h in findall(H.scenarioId .== s)) <= C_j[j] * y[j]) =#
                                    
        @constraint(mip, c5[j=eachindex(J), s=eachindex(S), t=eachindex(T)],
        L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h]
                            for h in sc_fps_range[s]) <= C_j[j] * y[j])
        #constraint (6)
        @constraint(mip, c6[j=eachindex(J), s=eachindex(S), t=2:length(T)],
                        L_ts[j, t, s] 
                        == 
                        L_ts[j, t-1, s] 
                        +
                        sum((λ[h, j, t] - b[h, j, t-1]) * u_h[h] for h in sc_fps_range[s]))

        #constraint (7)
        @constraint(mip, c7[j=eachindex(J), s=eachindex(scenarios)], L_ts[j, 1, s] == L_0[j])

        #constraint (8)
        @constraint(mip, c8_1[j=eachindex(J), t=eachindex(T), s=eachindex(scenarios)], L_ts[j, t, s] <= C_j[j] * y[j])
        #@constraint(mip, c8_2[j=eachindex(J), t=eachindex(T), s=eachindex(scenarios)], L_ts[j, t, s] >= 0)

        #constraint (9)
        #@constraint(mip, c9[j=eachindex(J)], L_0[j] >= 0)

        #constraint (10) and (11) are precised whene defining the variables


        # the objective function
        @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in findall(H.scenarioId .== s)) for s in eachindex(S)) -
                                sum(f_j[j] * y[j] for j in eachindex(J)) / cost_factor -
                                g * sum(L_0[j] for j in eachindex(J)) / cost_factor)


        #save the programme
        write_to_file(mip, mip_file_path)
    end
    set_optimizer(mip, Gurobi.Optimizer)
    set_time_limit_sec(mip, 3600.0)
    set_silent(mip)

    @info "solving the mip ..."
    optimize!(mip)

    ter_stat = termination_status(mip)
    if ter_stat in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.TIME_LIMIT]
        if ter_stat == MOI.TIME_LIMIT
            @warn "Time limit reached"
        end
        obj_val = objective_value(mip)
        cpu_time = solve_time(mip)
        sol_df = DataFrame(x=name.(all_variables(mip)), val=JuMP.value.(all_variables(mip)))
    
        sol = create_solution_for_simulation(sol_df, H)
        
        # return
        return obj_val, sol, cpu_time
    end
   
    return Inf, nothing, nothing
    
end
"""
    create_MIP(scenarios::Vector{Scenario}; save_mip_file = true, genarated_scenarios = false,
         mip_file_path="", costs_factors_list = [], mip_file_root = "")

Create the MIP model for the given scenarios.
If `save_mip_file` is true, the MIP model will be saved in the file `mip_file_path` 
    (or `mip_file_root`_cf.mof.json if 'costs_factors_list' is not empty) . 
If `genarated_scenarios` is true, the file where to save the mip will contain '_generated_' suffix. 
If `mip_file_root` is not empty, the costs_factors_list will mandatory be not empty and for each cost factor in the list, 
the objective function will be changed and the resulted mip will be saved in 'mip_file_root'_\$(cf)_cf.mof.json
"""
function create_MIP(scenarios::Vector{Scenario}; save_mip_file = true, genarated_scenarios = false,
         mip_file_path="", costs_factors_list = [], mip_file_root = "")
    global H
    global q_s
    global sc_fps_range
    # we make a convention that the file is names as "E_carsharing_mip_1_to_$(nbr_of_scenarios).mof.json"
    # so is we wanted to solve 5 scenarions so the name will be  "E_carsharing_mip_1_to_5.mof.json"
    if mip_file_path == ""
        mip_file_suffix = genarated_scenarios ? "_generated_$(length(scenarios))" : "1_to_$(length(scenarios))"
        mip_file_path = "Data/MIP/programs_file/E_carsharing_mip_$(mip_file_suffix)_$(maximum_walking_time)_walking_time.mof.json"
    end
   
    # The mixed integer program 
    mip = Model()

    # we need H variable in both case so we define it here
    # add column to the different scenarios to store the scenario id
    for i in eachindex(scenarios)
        scenarios[i].request_list.scenarioId = i .* ones(Int, nrow(scenarios[i].request_list)) 
        scenarios[i].feasible_paths.scenarioId = i .* ones(Int, nrow(scenarios[i].feasible_paths)) 
    end
    
    #Hs = [scenarios[i].feasible_paths  for i in eachindex(scenarios)]# all feasible paths
    #H = vcat(Hs...)
    H = vcat([scenarios[i].feasible_paths  for i in eachindex(scenarios)]...)
    
    # feasible paths range (to avoid filtring at each time)
    sc_fps_range= Array{UnitRange, 1}(undef, length(scenarios))
    curent_index = 1
    for sc in eachindex(scenarios)
        #get the range of feasible paths for the corresponding scenario
        sc_fps_range[sc] = UnitRange(curent_index, curent_index + nrow(scenarios[sc].feasible_paths) - 1)
        curent_index += nrow(scenarios[sc].feasible_paths)
    end
        
    @info "[Create MIP]: Start ..."
    #set the different variables for the MIP and try to preserve the same names as in the paper
    #K_s = [scenarios[i].request_list for i in eachindex(scenarios)]
    #K = vcat(K_s...) #set of all request in all scenarios
    K = vcat([scenarios[i].request_list for i in eachindex(scenarios)]...)
    #Δ_k = [get_trip_duration(req.ON, req.DN) / time_slot_length * 1.1 for req in eachrow(K)] #time threshold of the trip duration for each request
    
    S = collect(1:length(scenarios)) 
    # at this stage we are not taking in consideration the probability of each scenario
    global q_s = 1 / length(scenarios) .* ones(length(scenarios)) 

    @info "[Create MIP]: set the variables ..."
    #generate the parameters
    b = SparseArray{Bool}(undef, nrow(H), length(J), length(T)) 
    μ = SparseArray{Bool}(undef, nrow(H), length(J), length(T)) 
    λ = SparseArray{Bool}(undef, nrow(H), length(J), length(T))
    
    for h in 1:nrow(H)
        # set bata
        start_station_id = stations_idx[H.origin_station[h]]
        destination_station_id = stations_idx[H.destination_station[h]]
        starting_time_slot = H.start_driving_time[h]
        b[h, start_station_id, starting_time_slot] = 1
        #@benchmark b[h, stations_idx[H.origin_station[h]], H.start_driving_time[h]] = 1
        
        # set μ
        start_charging_time_slot = H.arriving_time[h]
        end_charging_time_slot = H.arriving_time[h] + ceil(Int, 
                            (δ_ij[start_station_id, destination_station_id] / charging_rate) / 60 / time_slot_length)
        μ[h, destination_station_id, start_charging_time_slot:end_charging_time_slot] .= 1
        
        # set λ
        λ[h, destination_station_id, end_charging_time_slot+1] = 1
    end
    
    # define the variables
    @variable(mip, u_h[1:nrow(H)], Bin) # u_h 1 if the trip in H is used 0 otherwise
    @variable(mip, L_ts[eachindex(J), eachindex(T), eachindex(S)] >= 0, Int)
    @variable(mip, L_0[eachindex(J)] >=0 , Int)
    @variable(mip, y[eachindex(J)], Bin)
    
    # List constraints
    @info "[Create MIP]: set the constraints ..."
    #constraints (2)
    #= @constraint(mip, c2[i=1:nrow(K)], 
                sum(u_h[j] for j in findall((H.req .== K.reqId[i]) 
                .& 
                (H.scenarioId .== K.scenarioId[i]))) <= 1) =#

    @constraint(mip, c2[i=1:nrow(K)], sum(u_h[j] for j in (K.fp[i] .+ (1000 * (K.scenarioId[i] - 1 ))) ) <= 1)

    #constraints (3)
    @constraint(mip, c33_1[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.origin_station[i]]])
    @constraint(mip, c33_2[i=1:nrow(H)], u_h[i] <= y[stations_idx[H.destination_station[i]]])
    
    #constraint (4)
    #=  
    @constraint(mip, c4[j=eachindex(J), s=eachindex(S), t=eachindex(T)], 
                    sum(b[h, j, t] * u_h[h] for h in sc_fps_range[s]) <= L_ts[j, t, s])
    =#
    
    @info "[Create MIP]: creating constraints 4"
    for j in eachindex(J) 
        for s in eachindex(S) 
            for t in eachindex(T)
                #non zero values of b
                bnz = [n[1] for n in nonzero_keys(b[sc_fps_range[s], j, t])]
                length(bnz) > 0 && @constraint(mip, [[j], [s], [t]], sum(b[h, j, t] * u_h[h] for h in bnz) <= L_ts[j, t, s])
                #@constraint(mip, [[j], [s], [t]], sum(b[h, j, t] * u_h[h] for h in sc_fps_range[s]) <= L_ts[j, t, s])

            end
        end
    end
          
    #constraint (5)
    #= 
    @constraint(mip, c5[j=eachindex(J), s=eachindex(S), t=eachindex(T)],
    L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h]
                        for h in sc_fps_range[s]) <= C_j[j] * y[j])
    =#
    
    @info "[Create MIP]: creating constraints 5"
    for j in eachindex(J)
        for s in eachindex(S)
            for t in eachindex(T)
                # get the indices of the non zero values of μ - b
                μbnz = [n[1] for n in union(nonzero_keys(b[sc_fps_range[s], j, t]), nonzero_keys(μ[sc_fps_range[s], j, t]))]
                length(μbnz) > 0 &&@constraint(mip, [[j], [s], [t]], L_ts[j, t, s] + sum((μ[h, j, t] - b[h, j, t]) * u_h[h]
                            for h in μbnz) <= C_j[j] * y[j])
                                                    
            end
        end
    end
    #constraint (6)
    #= @constraint(mip, c6[j=eachindex(J), s=eachindex(S), t=2:length(T)],
                    L_ts[j, t, s] 
                    == 
                    L_ts[j, t-1, s] 
                    +
                    sum((λ[h, j, t] - b[h, j, t-1]) * u_h[h] for h in sc_fps_range[s])) 
    =#
    @info "[Create MIP]: creating constraints 6"
    for j in eachindex(J)
        for s in eachindex(S)
            for t in 2:length(T)
                # get the indices of the non zero values of μ - b
                μbnz = [n[1] for n in union(nonzero_keys(b[sc_fps_range[s], j, t]), nonzero_keys(μ[sc_fps_range[s], j, t]))]
                
                @constraint(mip, [[j], [s], [t]], L_ts[j, t, s] == L_ts[j, t-1, s] 
                            +
                            sum((λ[h, j, t] - b[h, j, t-1]) * u_h[h] for h in μbnz))
            end
        end
    end

    # constraint (7)
    # @constraint(mip, c7[j=eachindex(J), s=eachindex(scenarios)], L_ts[j, 1, s] == L_0[j])
    @info "[Create MIP]: creating constraints 7"
    for j in eachindex(J)
        for s in eachindex(S)
            @constraint(mip, [[j], [s]], L_ts[j, 1, s] == L_0[j])
        end
    end
    #constraint (8)
    #@constraint(mip, c8_1[j=eachindex(J), t=eachindex(T), s=eachindex(scenarios)], L_ts[j, t, s] <= C_j[j] * y[j])
    @info "[Create MIP]: creating constraints 8"
    for j in eachindex(J)
        for s in eachindex(S)
            for t in eachindex(T)
                @constraint(mip, [[j], [t], [s]], L_ts[j, t, s] <= C_j[j] * y[j])
            end
        end
    end

    #@constraint(mip, c8_2[j=eachindex(J), t=eachindex(T), s=eachindex(scenarios)], L_ts[j, t, s] >= 0)
    # did it whene defining the variables
    
    #constraint (9)
    #@constraint(mip, c9[j=eachindex(J)], L_0[j] >= 0) did it in the definition of the variables

    #constraint (10) and (11) are precised whene defining the variables
    
    # the objective function
    @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in sc_fps_range[s]) for s in eachindex(S)) -
                            sum(f_j[j] * y[j] for j in eachindex(J)) / cost_factor -
                            g * sum(L_0[j] for j in eachindex(J)) / cost_factor)


    #check if costs factors list is not empty so we create an objective function for each cost factor
    if isempty(costs_factors_list)
        #save the programme
        save_mip_file && write_to_file(mip, mip_file_path)
    else
        for cf in costs_factors_list
            # set the new objective
            @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in sc_fps_range[s]) for s in eachindex(S)) -
                            sum(f_j[j] * y[j] for j in eachindex(J)) / cf -
                            g * sum(L_0[j] for j in eachindex(J)) / cf)
            
            # the file path 
            mip_file_path = string(mip_file_root,"_$(cf)_cf.mof.json")
            #save the programme
            save_mip_file && write_to_file(mip, mip_file_path)
        end
    end
end

""" f
    Create_MIPs_for_Generated_scenarios()

Create the MIPs for the generated scenarios. for different walking time and number of requests, and costs factors
"""
function create_MIPs_for_Generated_scenarios()
    if !work_with_time_slot
        work_with_time_slot = true
        @warn "The working with time slot is set to true!"
    end

    #first the generated Scenarios
    #list of parameters
    nbr_requests_list = [1000, 2000, 3000, 5000, 10000]
    walking_time_list = [5, 6, 7, 8, 10, 15]
    costs_factors_list = [10^4, 10^5, 10^6]
    
    #path where generated scenarios are saved
    generated_scs_folder_path = project_path("Data/generated_scenario/scenario_txt_files")

    # Dict to save time for creating the MIP
    time_df = DataFrame(params =Tuple[], mip_creation_time = Float64[])
    Threads.@threads for (nr, wt) in Iterators.product(nbr_requests_list, walking_time_list)
        @info "[Create MIP]: number of request $(nr) and walking time $(wt)"
        
        file_path = string(generated_scs_folder_path, "/scenario_", nr,"_requests.txt")

        global maximum_walking_time = wt

        #initialize the scenario
        scenario = initialize_scenario(file_path)

        #solve the MIP
        TT = @elapsed create_MIP([scenario], save_mip_file = true, costs_factors_list = costs_factors_list,
             mip_file_root = "Data/MIP/programs_file/E_carsharing_mip_generated_$(nr)_requests_$(wt)_walking_time")
        #= for cf in costs_factors_list
            @info "[Create MIP]: set objective for cost factor $(cf)"
            @objective(mip, Max, sum(q_s[s] * sum(H.Rev[h] * u_h[h] for h in sc_fps_range[s]) for s in eachindex(S)) -
                            sum(f_j[j] * y[j] for j in eachindex(J)) / cf -
                            g * sum(L_0[j] for j in eachindex(J)) / cf)
            
            # the mip file path where the MIP model will be saved
            mip_file_path = "Data/MIP/programs_file/E_carsharing_mip_generated_$(nr)_requests_$(wt)_walking_time_$(cf)_CF_.mof.json"
            #save the programme
            write_to_file(mip, mip_file_path)

            #add TT to the time dict
            time_dict[(nr, wt, cf)] = TT
        end =#
        for cf in costs_factors_list 
            push!(time_df, [(nr, wt, cf), TT])
        end
    end
    CSV.write("Data/MIP/programs_file/MIP_generated_scenarios_creation_time.csv", time_df)

end
"""
    create_MIPs_for_Ci()

for each Ci, create a MIP with respect to the parameters 
"""
function create_MIPs_for_Ci()
    # parameters for the experiments
    scenarios_sets = ["C1", "C2", "C3", "C4"]
    nbr_scenario_list = [100, 200]
    nbr_requests_list = [1000]
    walking_time_list = [5, 6, 7, 8, 9, 10, 15]
    costs_factors = [10^5, 10^6]

    time_df = DataFrame(params =Tuple[], mip_creation_time = Float64[])
    
    Threads.@threads for (set, ns, nr, wt) in Iterators.product(scenarios_sets, nbr_scenario_list, nbr_requests_list, walking_time_list)
        @info "[create_MIPs_for_Ci]: set = $set, number of scenarios = $ns, number of requests = $nr, walking time = $wt, cost_factor = $cf"
        
        #set the  global variables 
        global maximum_walking_time = wt
        global cost_factor = cf

        #prepare the scenarios
        scenarios = [initialize_scenario(project_path("Data/Instances/$set/scenario_txt_files/Output1000_$(set)_$(i).txt")) for i in 1:ns]
        
        # create the MIP and get the time of creation
        TT = @elapsed create_MIP(scenarios, save_mip_file = true, costs_factors_list = costs_factors_list,
            mip_file_root = "Data/MIP/programs_file/E_carsharing_mip_$(set)_$(ns)_scenarios_$(nr)_requests_$(wt)_walking_time")
  
        for cf in costs_factors_list 
            push!(time_df, [(nr, wt, cf), TT])
        end
        
    end
    CSV.write("Data/MIP/programs_file/MIP_creation_time.csv", time_df)
end
