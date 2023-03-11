
global results_folder = project_path("results")
#= 
    this function has the role of validating the simulation model by comparing the results of the simulation
    with the results of the mixed Integer programming model.
    - for a set of scenarios we create a mixed integer programming model and solve it using Gurobi.
    - then we got the solutions and evalute them using the simulation
    - finaly we compare between the results of the simulation and the results of the mixed integer programming model. 
=#
function validate_simulation_model()
    #we consider onlt ten scenarios
    scenario_ids = collect(1:10)
    
    #define different results variables
    mip_objective_values = Float64[]
    sim_objective_values = Float64[]
    solutions = Solution[]
    
    are_equal = Bool[]
    difference = Float64[]
    #run get_solutions
    for sc_id in scenario_ids
        @info " working with scenario N° $sc_id"
        #initiualize the scenario 
        initialize_scenarios([sc_id])
        sc = scenario_list[1]

        #solve the MIP
        mip_obj, sol = solve_using_mixed_integer_program([sc], mip_file_path = "Data/MIP/programs_file/E_carsharing_mip_$(sc_id).mof.json")
        
        #evaluate the solution using the simulation
        sim_obj = E_carsharing_sim(sol)

        #put the results in the corresponding variables
        push!(mip_objective_values, mip_obj)
        push!(sim_objective_values, sim_obj)
        push!(are_equal, isapprox(mip_obj, -1 * sim_obj, atol=1e-6))
        push!(difference, abs(mip_obj + sim_obj))
        
        #show a warning if the results are not equal
        if difference[sc_id] != 0.0
            message = string("the results are not equal for scenario $sc_id", are_equal[sc_id] ? ". However, they are almost equal" : "") 
            @warn message
        end
        push!(solutions, sol)
    end

    #save the solutions and objective values
    result_folder_for_this_experiment = string(results_folder, "/validation")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)
    objective_values_file_path = string(result_folder_for_this_experiment, "/objective_values.jls")
    solutions_file_path = string(result_folder_for_this_experiment, "/solutions.jls")

    #save to files
    serialize(objective_values_file_path, [mip_objective_values, sim_objective_values])
    serialize(solutions_file_path, solutions)

end

#= 
    this is reproduction for the first experiment the article [1].
    the aim is to have an overview of the pre-processing function how many requests are accesible 
    for different maximum walking times, how many feasible paths are generated and the CPU time 
    at each combination of number of requests and walking time (β_w).
    The result for this experiment in the article [1] is presented in TABLE I.

    references:
    [2]: Hatice Çalık, Bernard Fortz. A Benders decomposition method for locating stations in a one-way electric 
         car sharing system under demand uncertainty, Transportation Research Part B: Methodological, Volume 125,
         2019, Pages 121-150, ISSN 0191-2615
=#
function preprocessing_experiment()
    # global variables
    global maximum_walking_time
    global all_requests_list

    if work_with_time_slot
        work_with_time_slot = false
        @warn "The working with time slot is set to false"
    end

    # the path for the results folder for this experiment (created if doesn't exist)
    # the name preprocessing_2017 is because this experiment is the reproduction of the experiment on the article of 2017
    result_folder_for_this_experiment = string(results_folder, "/preprocessing_2017")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment,"/PS_", now(), ".csv")

    all_station = get_potential_locations()

    nbr_requests_list = [1000 , 2000, 3000, 5000, 10000]
    walking_time_list = [5, 6, 7, 8, 10, 15]
    
    #the results will be saved in a dataframe
    results_as_df = DataFrame(K = Int64[], β_w = Int64[], K_a = Int64[], H = Int64[], PP_time =[])
    
    for (nbr_requests, wt) in Iterators.product(nbr_requests_list, walking_time_list)
        # we are going to use the generated scenarios
        # set the file path
        curr_sc_path = "Data/generated_scenario/scenario_$(nbr_requests)_requests.txt"
        
        @info "($nbr_requests,$wt)  is being tested ..."
        curr_requests_list = requests_as_dataframe(curr_sc_path)
        
        #run the preprocessing_function
        curr_pp_time = @elapsed afp = get_feasible_paths(curr_requests_list, all_station, wt)
        #results traitement
        # 1-the accepted requets (|K_a|)
        curr_k_a =  length(unique(afp.req))
        # 2- feasible paths size (|H|)
        curr_H = nrow(afp)
        # push the results
        push!(results_as_df, [nbr_requests, wt, curr_k_a, curr_H, curr_pp_time])        
    end
    CSV.write(results_save_path, results_as_df)
end

#= 
    As the previous experiment, we seek to see the number of accessible requests 
    as a function of Number of scenarios and maximum walking time.
    this experiment is the reproduction of the experiment conducted in [2] and results are presented in Fig 2 page 127
    references:
    [2]: Hatice çalik, Bernard fortz. Location of Stations in One-Way Electric Car sharing System.
         IEEE Symposium on Cpmputer and Communications, 2017 Heraklion, Greec. hal-01665609
 =#
function preprocessing_experiment2019()
    # global variables
    global all_requests_list
    
    # as we are interested by only the feasibles paths, we don't need to work with time slot
    if work_with_time_slot
        work_with_time_slot = false
        @warn "The working with time slot is set to true!"
    end
    # the path for the results folder for this experiment (created if doesn't exist)
    # the name preprocessing_2019 is because this experiment is the reproduction of the experiment on the article of 2019
    result_folder_for_this_experiment = string(results_folder, "/preprocessing_2019")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment,"/PS_scenarios_", now(), ".csv")

    all_station = get_potential_locations()

    nbr_scenarions_list =[100, 200]
    walking_time_list = [5, 6, 7, 8, 10, 15]
    results_as_df = DataFrame(K = Int64[], β_w = Int64[], K_a = Int64[], H = Int64[], PP_time =[])

    #initialize scenarios to get the requests list
    initialize_scenarios(collect(1:maximum(nbr_scenarions_list)))
    for (scenario_number, wt) in Iterators.product(nbr_scenarions_list, walking_time_list)
        @info "($scenario_number,$wt)  is being tested ..."
        curr_requests_list = vcat([sc.request_list for sc in scenario_list[1:scenario_number]]...)
        
        #run the preprocessing_function
        curr_pp_time = @elapsed afp = get_feasible_paths(curr_requests_list, all_station, wt)
        #results traitement
        # 1-the accepted requets (|K_a|)
        curr_k_a =  length(unique(afp.req))
        # 2- feasible paths size (|H|)
        curr_H = nrow(afp)
        # push the results
        push!(results_as_df, [scenario_number, wt, curr_k_a, curr_H, curr_pp_time])        
    end
    CSV.write(results_save_path, results_as_df)
end


#= construct the scenarios for the next experiment 
    basically make  scenarios with n requests by concatinating the different request from the existing requests of 1000
=#
function construct_scenario_with_different_size(sizes::Vector{Int64})
    #the path where the file will be stored
    folder_path = project_path("Data/generated_scenario")
    # count the number of scenarions with 1000 that we have to initialize
    total_number_of_requests = maximum(sizes)
    nbr_files_to_read = ceil(Int64, total_number_of_requests/1000)
    
    #read the files
    requests_ids_list = vcat([readlines(scenarios_paths[i]) for i in 1:nbr_files_to_read]...)

    for s in sizes
        file_path = string(folder_path, "/scenario_", s,"_requests.txt")
        writedlm(file_path, requests_ids_list[1:s])
    end
end
#=
    this function solve the Mixed Integer programms and return the best value found
    the aim is to compare afterwords the result with the Hyper heuristic frame work 
    N.P: another same experiment will be done in Java counterpart to compare the results
    The results can not be compared to TABLE II - IV in article [1] as we don't know what are the requests used. Nevertheless 
    we can have an idea about the values ( it is always benificial to see if we are in the same range of values)
 =#
 
function mixed_integer_programming_experiment()
    global work_with_time_slot
    # we will try to use the same variables name as the paper 
    if !work_with_time_slot
        work_with_time_slot = true
        @warn "The working with time slot is set to true!"
    end

    generated_scs_folder_path = project_path("Data/generated_scenario")
    result_folder_for_this_experiment = string(results_folder, "/mixed_integer_programming_experiment")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment,"/MIP_", now(), ".csv")
    
    #list of parameters
    nbr_requests_list = [1000, 2000, 3000, 5000, 10000]
    walking_time_list = [5, 6, 7, 8, 10, 15]
    costs_factors_list = [10^4, 10^5, 10^6]

    results_as_df = DataFrame(CF=[], K=[], β_w=[], PF_Opt=[], J_bar=[], K_bar=[], solver_time=[], total_time =[])
    for (cf, nr, wt) in Iterators.product(costs_factors_list, nbr_requests_list, walking_time_list)
        @info "[MIP experiment Cost factor = $cf]: solving a MIP with $nr requests and walking time equal to $wt..."
        
        file_path = string(generated_scs_folder_path, "/scenario_", nr,"_requests.txt")
       
        global maximum_walking_time = wt
        global cost_factor = cf

        #initialize the scenario
        scenario = initialize_scenario(file_path, check_file=false)

        # the mip file path where the MIP model will be saved
        mip_file_path = "Data/MIP/programs_file/E_carsharing_mip_generated_$(nr)_requests_$(wt)_walking_time_$(cf)_CF_.mof.json"
        
        #solve the MIP
        TT = @elapsed obj, sol, solve_time = solve_with_mixed_integer_program([scenario], mip_file_path)
        if obj == -Inf
            push!(results_as_df, [cf, nr, wt, obj, missing, missing, solve_time, TT])
        else
            push!(results_as_df, [cf, nr, wt, obj, sum(sol.open_stations_state), sum(sol.selected_paths[1]), solve_time, TT])

            #save the sol file
            save_sol_path = string(generated_scs_folder_path, "/serialized_solutions/sol_", nr,"_requests.jls")
            serialize(save_sol_path, sol)
        end
    end
    CSV.write(results_save_path, results_as_df)
end

#create a function to generate feasible paths for the generated scenarios
function generate_feasible_paths_for_generated_scenarios()
    #list of walking times
    walking_time_list = [5, 6, 7, 8, 10, 15]
    all_station = potential_locations
    gen_sce_paths = filter!(x -> startswith(x, "scen"), readdir("Data/generated_scenario"))
    #loop over scenarios_paths
    for path in gen_sce_paths 
        requests = requests_as_dataframe("Data/generated_scenario/$(path)")
        
        for wt in  walking_time_list
            @info "generating feasible paths for $(path) with walking time equal to $wt"
            #get the feasible paths
            fps = get_feasible_paths(requests, all_station, wt)
            
            #save the feasible paths as CSV file
            save_path = string("Data/feasible_paths/paths_generated_", path,"_", wt, "min_walking_time.csv")

            #save the file
            CSV.write(save_path, fps)

        end
    end
end

