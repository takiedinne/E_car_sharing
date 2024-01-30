
global results_folder = project_path("results")
""" 
    this function has the role of validating the simulation model by comparing the results of the simulation
    with the results of the mixed Integer programming model.
    - for a set of scenarios we create a mixed integer programming model and solve it using Gurobi.
    - then we got the solutions and evalute them using the simulation
    - finaly we compare between the results of the simulation and the results of the mixed integer programming model. 
"""
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
        mip_obj, sol = solve_using_mixed_integer_program([sc], mip_file_path="Data/MIP/programs_file/E_carsharing_mip_$(sc_id).mof.json")

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

""" 
    this is reproduction for the first experiment the article [1].
    the aim is to have an overview of the pre-processing function how many requests are accesible 
    for different maximum walking times, how many feasible paths are generated and the CPU time 
    at each combination of number of requests and walking time (β_w).
    The result for this experiment in the article [1] is presented in TABLE I.

    references:
    [2]: Hatice Çalık, Bernard Fortz. A Benders decomposition method for locating stations in a one-way electric 
         car sharing system under demand uncertainty, Transportation Research Part B: Methodological, Volume 125,
         2019, Pages 121-150, ISSN 0191-2615
"""
function preprocessing_experiment()
    println("[preprocessing_experiment 2017]: start ...")
    # global variables
    global maximum_walking_time
    global all_requests_list
    global work_with_time_slot

    if work_with_time_slot
        work_with_time_slot = false
        @warn "The working with time slot is set to false"
    end

    # the path for the results folder for this experiment (created if doesn't exist)
    # the name preprocessing_2017 is because this experiment is the reproduction of the experiment on the article of 2017
    result_folder_for_this_experiment = string(results_folder, "/preprocessing_2017")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file

    #now_as_str = Dates.format(now(), "yyyy-mm-ddTHH_MM_SS")
    results_save_path = project_path(string(result_folder_for_this_experiment, "/PS_", now_as_str, ".csv"))

    all_station = get_potential_locations()

    nbr_requests_list = [1000, 2000, 3000, 5000, 10000]
    walking_time_list = [5, 6, 7, 8, 10, 15]

    #the results will be saved in a dataframe
    results_as_df = DataFrame(K=Int64[], β_w=Int64[], K_a=Int64[], H=Int64[], PP_time=[])

    for (nbr_requests, wt) in Iterators.product(nbr_requests_list, walking_time_list)
        # we are going to use the generated scenarios
        # set the file path
        curr_sc_path = project_path("Data/generated_scenario/scenario_txt_files/scenario_$(nbr_requests)_requests.txt")

        @info "[preprocessing experiment 2017]: ($nbr_requests,$wt)  is being tested ..."
        curr_requests_list = requests_as_dataframe(curr_sc_path)

        #run the preprocessing_function
        curr_pp_time = @elapsed afp = get_feasible_paths(curr_requests_list, all_station, wt)
        #results traitement
        # 1-the accepted requets (|K_a|)
        curr_k_a = length(unique(afp.req))
        # 2- feasible paths size (|H|)
        curr_H = nrow(afp)
        # push the results
        push!(results_as_df, [nbr_requests, wt, curr_k_a, curr_H, curr_pp_time])
    end
    CSV.write(results_save_path, results_as_df)
    @info "[preprocessing Exp 2017]: The experiment is finished !"
end


"""
construct the scenarios for the next experiment 
    basically make  scenarios with n requests by concatinating the different request from the existing requests of 1000
"""
function construct_scenario_with_different_size(sizes::Vector{Int64})
    #the path where the file will be stored
    folder_path = project_path("Data/generated_scenario")
    # count the number of scenarions with 1000 that we have to initialize
    total_number_of_requests = maximum(sizes)
    nbr_files_to_read = ceil(Int64, total_number_of_requests / 1000)

    #read the files
    requests_ids_list = vcat([readlines(scenarios_paths[i]) for i in 1:nbr_files_to_read]...)

    for s in sizes
        file_path = string(folder_path, "/scenario_", s, "_requests.txt")
        writedlm(file_path, requests_ids_list[1:s])
    end
end
"""
    this function solve the Mixed Integer programms and return the best value found
    the aim is to compare afterwords the result with the Hyper heuristic frame work 
    N.P: another same experiment will be done in Java counterpart to compare the results
    The results can not be compared to TABLE II - IV in article [1] as we don't know what are the requests used. 
    Nevertheless, we can have an idea about the values ( it is always benificial to see if we are in the same range 
    of values)
"""
function solve_generated_scenario_using_Gurobi()
    global work_with_time_slot
    # we will try to use the same variables name as the paper 
    if !work_with_time_slot
        work_with_time_slot = true
        @warn "The working with time slot is set to true!"
    end

    generated_scs_folder_path = project_path("Data/generated_scenario")
    result_folder_for_this_experiment = string(results_folder, "/GUROBI_for_generated_scenarios_2017")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/MIP_", now(), ".csv")

    #list of parameters
    nbr_requests_list = [1000] #= , 2000, 3000, 5000, 10000 =#
    walking_time_list = [5] #= , 6, 7, 8, 10, 15 =#
    costs_factors_list = [10^5] #= , 10^6 =#

    # we store theresults in a dataframe object
    results_as_df = DataFrame(CF=[], K=[], β_w=[], PF_Opt=[], J_bar=[], K_bar=[], solver_time=[], total_time=[], termination_status=[])
    for (cf, nr, wt) in Iterators.product(costs_factors_list, nbr_requests_list, walking_time_list)
        @info "[GUROBI Experiment]: cf = $cf, nr = $nr, wt = $wt is being solved ..."

        #noramally all the MIPs are generated befor. However, we check again and if they don't exist we generate them
        # the mip file path where the MIP model will be saved
        mip_file_path = project_path("Data/MIP/programs_file/E_carsharing_mip_generated_$(nr)_requests_$(wt)_walking_time_$(cf)_cf.mof.json")
        # path for the current Scenario
        curr_sc_path = project_path("Data/generated_scenario/scenario_txt_files/scenario_$(nr)_requests.txt")
        #initialize the scenario
        scenario = initialize_scenario(curr_sc_path)
        if !isfile(mip_file_path)
            @warn "The MIP file is not found, it will be generated ..."
            file_path = string(generated_scs_folder_path, "/scenario_txt_files/scenario_", nr, "_requests.txt")

            # set the global variables
            global maximum_walking_time = wt
            global cost_factor = cf

            #create the MIPs and save them
            create_MIP([scenario], save_mip_file=true, costs_factors_list=costs_factors_list,
                mip_file_root="Data/MIP/programs_file/E_carsharing_mip_generated_$(nr)_requests_$(wt)_walking_time")
        end
        #=  
            #here we are 100% sure that the MIP file exists. solve_using_mixed_integer_program as well can creat MIP
            # Nevertheless, I decided to creat them manually as foreach pair of nr, wt we can create 3 MIPs (one for each cf) 
            # by changing only the Objective function 
        =#

        #solve the MIP 
        TT = @elapsed obj, sol, solve_duration, ter_stat = solve_using_mixed_integer_program([scenario], mip_file_path=mip_file_path)
        if obj == Inf
            push!(results_as_df, [cf, nr, wt, obj, missing, missing, solve_time, TT, ter_stat])
        else
            push!(results_as_df, [cf, nr, wt, obj, sum(sol.open_stations_state), sum(sol.selected_paths[1]), solve_time, TT, ter_stat])

            #save the sol file
            save_sol_path = string(generated_scs_folder_path, "/serialized_solutions/sol_$(nr)_requests_$(wt)_walking_time_$(cf)_cf.jls")
            serialize(save_sol_path, sol)
        end
    end
    CSV.write(results_save_path, results_as_df)
end

"""
    create a function to generate feasible paths for the generated scenarios
"""
function generate_feasible_paths_for_generated_scenarios()
    #list of walking times
    walking_time_list = [5, 6, 7, 8, 10, 15]
    all_station = get_potential_locations()

    gen_sce_paths = filter!(x -> startswith(x, "scen"), readdir("Data/generated_scenario/scenario_txt_files"))
    #loop over scenarios_paths
    for path in gen_sce_paths
        requests = requests_as_dataframe("Data/generated_scenario/scenario_txt_files/$(path)")

        for wt in walking_time_list
            @info "generating feasible paths for $(path) with walking time equal to $wt"
            #get the feasible paths
            fps = get_feasible_paths(requests, all_station, wt)

            #save the feasible paths as CSV file
            save_path = string("Data/generated_scenario/feasible_paths/paths_generated_", path, "_", wt, "min_walking_time.csv")

            #save the file
            CSV.write(save_path, fps)

        end
    end
end

"""
    create the set of instancies C1, ..., C4
    the idea is to sample uniformly from 1:1000 a set of 200 scenarios, 4 times to 
    construct the scenarios sets C_i for i ∈ {1, ..., 4}
"""
function construct_scenarios_sets()
    #create the folder where the scenarios will be stored if not exists
    if !isdir("Data/Instances/C1")
        mkpath("Data/Instances/C1")
        mkpath("Data/Instances/C2")
        mkpath("Data/Instances/C3")
        mkpath("Data/Instances/C4")
    end

    #sample four time 200 scenarios from the 1000 scenarios
    for i in 1:4
        @info "constructing the set C$i"
        #sample 200 scenarios
        sampled_scenarios = sample(1:1000, 200, replace=false)

        #copy the corresponding files to the relative folder
        for j in eachindex(sampled_scenarios)
            sc = sampled_scenarios[j]
            src_path = "Data/Scenarios_1000_greaterthan2/Output1000_$sc.txt"

            #copy the file
            cp(src_path, "Data/Instances/C$i/Output1000_C$(i)_$(j).txt", force=true)
        end

    end
end

"""
    solve a set of scenarios  using Mixed Integer programming
    we try to reproduce the Fig3 in Article [2]
    references:
    [2]: Hatice Çalık, Bernard Fortz. A Benders decomposition method for locating stations in a one-way electric 
         car sharing system under demand uncertainty, Transportation Research Part B: Methodological, Volume 125,
         2019, Pages 121-150, ISSN 0191-2615
"""
function solve_Ci_set_with_MIP()
    #list of parameters
    nbr_scenario_list = [10, 100, 200]
    nbr_requests_list = [1000, 1000, 5000]
    walking_time_list = [5, 5, 5]
    costs_factors = [10^5, 10^5, 10^5]
    # make sure that we are working with time slot
    global work_with_time_slot
    if !work_with_time_slot
        work_with_time_slot = true
        @warn "The working with time slot is set to true!"
    end

    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/solve_ci_sets_with_MIP")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)
    !isdir(project_path("Data/MIP/solutions")) && mkpath(project_path("Data/MIP/solutions"))
    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/Ci_sets_MIP_Gurobi_", now(), ".csv")

    # parameters for the experiments    
    results_as_df = DataFrame(NS=Int64[], K=Int64[], β_w=[], PF_Opt=[], J_bar=Int64[], K_bar=Int64[], solver_time=[], total_time=[], terminal_status=[])
    CSV.write(results_save_path, results_as_df)
    #for (ns, nr, wt, cf) in Iterators.product(nbr_scenario_list, nbr_requests_list, walking_time_list, costs_factors)
    for i in eachindex(nbr_scenario_list)
        ns, nr, wt, cf = nbr_scenario_list[i], nbr_requests_list[i], walking_time_list[i], costs_factors[i]
        @info "[Ci MIP experiment 2019]: number of scenarios = $ns, number of requests = $nr, walking time = $wt, cost_factor = $cf"

        #set the  global variables 
        global maximum_walking_time = wt
        global number_of_requests_per_scenario = nr
        global cost_factor = cf

        # Mip file path
        mip_file_path = project_path("Data/MIP/programs_file/E_carsharing_mip_$(ns)_scenarios_$(nr)_requests_$(wt)_walking_time.mof.json")
        sol_file_path = project_path("Data/MIP/solutions/E_carsharing_mip_$(ns)_scenarios_$(nr)_requests_$(wt)_walking_time.jls")

        TT = @elapsed begin
            #prepare the scenarios
            #@info "inititialize scenarios ..."
            scenarios = [initialize_scenario(project_path("Data/Instances/C1_5000_500/scenario_txt_files/Output_$(i).txt"), nr) for i in 1:ns]

            # solve using MIP solver 
            #@info "solving the scenarios ..."
            obj, sol, cpu_time, solver_ter_state = solve_using_mixed_integer_program(scenarios, mip_file_path=mip_file_path)
        end
        #save the results
        if obj == Inf
            push!(results_as_df, [ns, nr, wt, obj, missing, missing, cpu_time, TT])
        else
            empty!(results_as_df)
            push!(results_as_df, [ns, nr, wt, obj, sum(sol.open_stations_state), sum([sum(sol.selected_paths[i]) for i in eachindex(scenarios)]), cpu_time, TT, solver_ter_state])

            CSV.write(results_save_path, results_as_df, append=true)
            #save the sol file
            serialize(sol_file_path, sol)
        end
    end
    #CSV.write(results_save_path, results_as_df)
    results_as_df
end


"""
    create a function to generate feasible paths for the generated scenarios
"""
function generate_feasible_paths_for_Ci()
    #list of walking times
    walking_time_list = [5, 10, 15] #= 6, 7, 8, 9, =#
    nbr_requests_list = [1000, 2000, 5000]

    @info "generating feasible paths for the set C1"
    #scenarios_paths
    scenarios_paths = project_path.("Data/Instances/C1_5000_500/scenario_txt_files/" .* filter!(x -> startswith(x, "Out"), readdir(project_path("Data/Instances/C1_5000_500/scenario_txt_files"))))

    #loop over scenarios_paths
    for (path, wt, nr) in Iterators.product(scenarios_paths, walking_time_list, nbr_requests_list)

        #set the maximum walking time
        global maximum_walking_time = wt
        global number_of_requests_per_scenario = nr

        #initialize the scenario (N.P: if it was initilized befor so the function will loadit directly)
        current_scenario = initialize_scenario(path)

        fps = current_scenario.feasible_paths

        #save feasible paths the file
        # Replace "scenario_txt_files" with "feasible_paths"
        save_fps_path = replace(path, "scenario_txt_files" => "feasible_paths")
        # Replace ".txt" with ".csv"
        save_fps_path = replace(save_fps_path, r"\.txt$" => "_$(number_of_requests_per_scenario)_requests_$(maximum_walking_time)_walking_time.csv")

        CSV.write(save_fps_path, fps)

    end

end

"""
    function check_the_influence_of_the_defferences_with_hatice_preprocessing()
    
    try to see the infulence of the deferences between the preprocessing of Taki and Hatice.
    for this experiment we have three type of speeds which affect the driving duration 
        {fixed_driving_speed, multiple_driving_speeds, time_dependent_driving_speeds}
    for the walking graph we have two types of graphs
        {Directed, undirected}
"""
function check_the_influence_of_the_defferences_with_hatice_preprocessing()
    #the parameters for the experiments
    speed_types = ["multiple_driving_speeds", "multiple_driving_speeds", "fixed_driving_speed", "fixed_driving_speed"] #= , "time_dependent_driving_speeds", "time_dependent_driving_speeds" =#
    walking_graph_types = ["Directed", "undirected", "Directed", "undirected"] #= , "Directed", "undirected" =#

    multiple_driving_speeds_list = [true, true, false, false] #= , false, false =#
    use_dynamic_speeds_list = [false, false, false, false] #= , true, true =#
    walking_graph_types_list = [MetaDiGraph, MetaGraph, MetaDiGraph, MetaGraph] #= , MetaDiGraph, MetaGraph =#


    #the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/preprocessing_Taki_VS_Hatice")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = "$(result_folder_for_this_experiment)/Taki_Vs_hatice.csv"

    # load the request of the scenario 
    scenario_requests = requests_as_dataframe(project_path("Data/Instances/C1_5000_500/scenario_txt_files/Output_1.txt"))

    # result object as dataframe
    results = DataFrame(speed_type=String[], walking_graph_type=String[], nbr_feasible_trips=Int64[], nbr_feasible_requests=Int64[])

    # try the get the feasible path for the scenario
    for i in eachindex(multiple_driving_speeds_list)

        #set the values for each parameters
        global multiple_driving_speeds = multiple_driving_speeds_list[i]
        global use_dynamic_speeds = use_dynamic_speeds_list[i]
        global graph_type = walking_graph_types_list[i]

        @info "[preprocessing Taki vs hatice]: multiple_driving_speeds = $multiple_driving_speeds, use_dynamic_speeds = $use_dynamic_speeds, graph_type = $graph_type"

        # load the adequate driving graph
        global manhaten_city_driving_graph = (multiple_driving_speeds) ? loadgraph(Manhatten_network_driving_graph_file, MGFormat()) : loadgraph(Manhatten_network_length_graph_file, MGFormat())
        global manhaten_city_length_graph = graph_type(loadgraph(Manhatten_network_length_graph_file, MGFormat()))

        #clean up the walking_time_list
        empty!(shortest_car_paths)
        empty!(shortest_walking_paths)

        #get the feasible paths
        feasible_paths = get_feasible_paths(scenario_requests, get_potential_locations(), 5)

        push!(results, [speed_types[i], walking_graph_types[i], nrow(feasible_paths), length(unique(feasible_paths.req))])

    end
    CSV.write(results_save_path, results)
end


""" 
    As the previous experiment, we seek to see the number of accessible requests 
    as a function of Number of scenarios and maximum walking time.
    this experiment is the reproduction of the experiment conducted in [2] and results are presented in Fig 2 page 127
    references:
    [2]: Hatice çalik, Bernard fortz. Location of Stations in One-Way Electric Car sharing System.
         IEEE Symposium on Cpmputer and Communications, 2017 Heraklion, Greec. hal-01665609
"""
function preprocessing_experiment2019()
    # Experimenents parameters
    nbr_scenarios_list = [100, 200, 300, 400, 500]
    walking_time_list = [5, 6, 7, 8, 9, 10, 15]
    scenario_size_list = [1000, 2000, 3000, 4000, 5000]

    global scenarios_paths = project_path.("Data/Instances/C1_5000_500/scenario_txt_files/" .* filter!(x -> startswith(x, "Out"), readdir(project_path("Data/Instances/C1_5000_500/scenario_txt_files"))))

    # global variables
    global all_requests_list
    global work_with_time_slot
    # as we are interested by only the feasibles paths, we don't need to work with time slot
    if work_with_time_slot
        work_with_time_slot = false
        @warn "The working with time slot is set to False"
    end
    # the path for the results folder for this experiment (created if doesn't exist)
    # the name preprocessing_2019 is because this experiment is the reproduction of the experiment on the article of 2019
    result_folder_for_this_experiment = string(results_folder, "/preprocessing_2019")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = project_path(string(result_folder_for_this_experiment, "/preprocessing_results", ".csv"))

    all_station = get_potential_locations()

    results_as_df = DataFrame(S=[], K=Int64[], β_w=Int64[], K_a=Int64[], H=Int64[], PP_time=[])

    @info "[preprocessing Exp 2019]: start ..."
    #initialize scenarios to get the requests list
    scenarios_requests_list = []
    global number_of_requests_per_scenario = maximum(scenario_size_list)
    for i in 1:maximum(nbr_scenarios_list)

        requests = requests_as_dataframe(scenarios_paths[i])
        requests.reqId = requests.reqId .+ (i - 1) * maximum(scenario_size_list)
        push!(scenarios_requests_list, requests)
    end

    for (scenario_number, wt, scenario_size) in Iterators.product(nbr_scenarios_list, walking_time_list, scenario_size_list)
        @info "[preperocessing Exp 2019]: |S| = $scenario_number, wt = $wt, nbr_requests = $scenario_size is being tested ..."

        curr_requests_list = vcat([scenarios_requests_list[i][1:scenario_size, :] for i in 1:scenario_number]...)

        #run the preprocessing_function
        curr_pp_time = @elapsed afp = get_feasible_paths(curr_requests_list, all_station, wt)
        #results traitement
        # 1-the accepted requets (|K_a|)
        curr_k_a = length(unique(afp.req))
        # 2- feasible paths size (|H|)
        curr_H = nrow(afp)
        # push the results
        push!(results_as_df, [scenario_number, scenario_size, wt, curr_k_a, curr_H, curr_pp_time])
    end
    CSV.write(results_save_path, results_as_df)
    @info "[preprocessing Exp 2019]: The experiment is finished !"
    results_as_df
end


function solve_single_scenario_using_gurobi()
    #list of parameters
    scenario_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    walking_time_list = [5, 10, 15]



    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/solve_single_scenarios_with_MIP")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)
    !isdir(project_path("Data/MIP/solutions")) && mkpath(project_path("Data/MIP/solutions"))

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/single_scenarios_MIP_Gurobi_", now(), ".csv")

    # parameters for the experiments    
    results_as_df = DataFrame(S_id=Int64[], β_w=[], PF_Opt=[], J_bar=Int64[], K_bar=Int64[], solver_time=[], total_time=[], terminal_status=[])

    for (sc_id, wt) in Iterators.product(scenario_list, walking_time_list)

        @info "[Ci MIP experiment 2019]: scenario N° $sc_id walking time = $wt"

        #set the  global variables 
        global maximum_walking_time = wt

        # Mip file path
        mip_file_path = project_path("Data/MIP/programs_file/E_carsharing_mip_scenario_$(sc_id)_requests_1000_walking_time_$(wt).mof.json")
        sol_file_path = project_path("Data/MIP/solutions/E_carsharing_mip_scenario_$(sc_id)_requests_1000_walking_time_$(wt).jls")

        TT = @elapsed begin
            #prepare the scenarios
            #@info "inititialize scenarios ..."
            scenarios = initialize_scenarios([sc_id])

            # solve using MIP solver 
            #@info "solving the scenarios ..."
            obj, sol, cpu_time, solver_ter_state = solve_using_mixed_integer_program(scenarios, mip_file_path=mip_file_path)
        end
        #save the results
        if obj == Inf
            push!(results_as_df, [ns, nr, wt, obj, missing, missing, cpu_time, TT])
        else
            push!(results_as_df, [sc_id, wt, obj, sum(sol.open_stations_state), sum([sum(sol.selected_paths[i]) for i in eachindex(scenarios)]), cpu_time, TT, solver_ter_state])

            #save the sol file
            serialize(sol_file_path, sol)
        end
    end
    CSV.write(results_save_path, results_as_df)
    results_as_df
end


function solve_single_scenario_using_SA()
    #global variables experiment related variables
    global rng
    trial_nbr = 10
    main_seed = 1905

    #variables related to simulated simulated_annealing
    T, T₀, I, α, β = 796.0, 20.0, 82, 0.98, 0.8

    #scenario parameters
    scenario_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    walking_time_list = [5]

    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/solve_single_scenario_with_SA")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/single_scenario_SA_", now(), ".csv")

    df = DataFrame(scenario_id=Int64[], trial_id=Int64[], sa_obj=[], cpu_time=[])
    CSV.write(results_save_path, df)# write the header

    for (sc_id, wt) in Iterators.product(scenario_list, walking_time_list)

        @info "[ECS using SA]: scenario N° $sc_id walking time = $wt"

        #set the  global variables 
        global maximum_walking_time = wt

        #initialize the scenario:
        initialize_scenarios([sc_id])
        #global request_feasible_trips_ids = [] #to get feasible trips for req i on scenario s : request_feasible_trips_ids[s][i]

        initial_solutions = [generate_random_solution() for _ in 1:trial_nbr]

        for i in 1:trial_nbr
            rng = MersenneTwister(main_seed + i)
            _, obj, sa_cpu = simulated_annealing(initial_solutions[i], T, T₀, α, I, β)
            empty!(df)
            push!(df, [sc_id, i, obj, sa_cpu])
            CSV.write(results_save_path, df, append=true)
        end

    end

end

function SA_params()
    global rng

    T_list = Float64[1061.95]
    T0_list = Float64[1.]
    α_list = [0.98] #= , 0.99, 0.998 =#
    I_list = [82]
    β_list = collect(0.5:0.01:0.9)

    main_seed = 1905

    trial_nbr = 1
    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/SA_params/Beta")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/SA_params", now(), ".csv")

    #parameters for the experiments    
    df = DataFrame(T=Float64[], T0=Float64[], α=Float64[], I=Float64[], β=Float64[], SA_Obj=Float64[], best_SA=Float64[], total_time=Float64[])
    CSV.write(results_save_path, df) #= , append = true =#
    #prepare the scenarios
    #@info "inititialize scenarios ..."
    initialize_scenarios([1])

    starting_sols = [generate_random_solution() for _ in 1:trial_nbr]
    for (T, T0, α, I, β) in Iterators.product(T_list, T0_list, α_list, I_list, β_list)
        @info "T = $T, T0 = $T0, α = $α, I = $I, β = $β"

        fit = 0.0
        best_fit = Inf
        TT = @elapsed for i in 1:trial_nbr
            rng = MersenneTwister(main_seed + i)
            _, obj, _ = simulated_annealing(starting_sols[i], T, T0, α, I, β)
            fit += obj
            obj < best_fit && (best_fit = obj)
        end

        empty!(df)
        push!(df, [T, T0, α, I, β, fit / trial_nbr, best_fit, TT / trial_nbr])
        CSV.write(results_save_path, df, append=true)
    end
end

function SA_get_T0()
   
    T = 1
    acceptance_rate = 0.8
    incrementing_factor = 1.1
    emperical_rate = 0.0 
    while emperical_rate < acceptance_rate
        @info "trying T = $T"
        T = T * incrementing_factor
        emperical_rate = mean([get_acceptance_rate(T) for _ in 1:10])
    end

    return T
end


function get_acceptance_rate(T)
    sc_id = 7
    initialize_scenarios([sc_id])
    
    sol = generate_random_solution()
    
    sol_fit = E_carsharing_sim(sol)

    empirical_acceptance_rate = 0.0
    total_number_of_neighbors = length(findall(sol.open_stations_state))

    #@info "trying T = $T"
    accepted_sol = 0

    # second neighbourhood
    open_stations = findall(sol.open_stations_state)

    for i in open_stations
        new_sol = deepcopy(sol)
        new_sol.open_stations_state[i] = false
        new_sol.initial_cars_number[i] = 0

        lost_requests = clean_up_trips!(new_sol, scenario_list, i)

        # step 4: reassign the lost requests if we can
        assigne_requests!(new_sol, scenario_list, lost_requests)

        new_fit = ECS_objective_function(new_sol)

        if new_fit < sol_fit
            accepted_sol += 1
        else
            if rand(rng) <= exp(-(new_fit - sol_fit) / T)
                accepted_sol += 1
            end
        end
    end

    #calculate the empirical_acceptance_rate
    empirical_acceptance_rate = accepted_sol / total_number_of_neighbors

    return empirical_acceptance_rate
end

function solve_single_scenario_using_GRA()
       
    #scenario parameters
    scenario_list_id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    walking_time_list = [5]

    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/solve_single_scenario_with_greedy_assign_requests")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/single_scenario_GAR_", now(), ".csv")

    df = DataFrame(scenario_id=Int64[], trial_id=Int64[], sa_obj=[], cpu_time=[])
    CSV.write(results_save_path, df)# write the header

    for (sc_id, wt) in Iterators.product(scenario_list_id, walking_time_list)
        
        @info "[ECS using SA]: scenario N° $sc_id walking time = $wt"

        #set the  global variables 
        global maximum_walking_time = wt

        #initialize the scenario:
        initialize_scenarios([sc_id])
        initialize_sim(Solution(), scenario_list[1])
        #global request_feasible_trips_ids = [] #to get feasible trips for req i on scenario s : request_feasible_trips_ids[s][i]

        _, obj, sa_cpu = greedy_assign_requests()
        empty!(df)
        push!(df, [sc_id, 1, obj, sa_cpu])
        CSV.write(results_save_path, df, append=true)
        
    end

end

function ruin_depth_exp()
    global ruin_depth
    global rng

    seed = 1905
 
    ruin_depth_list = collect(0.005:0.001:0.015)
    trial_nbr = 1
    scenario_list_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    # the folder where the results will be stored
    result_folder_for_this_experiment = string(results_folder, "/ruin_depth")
    !isdir(result_folder_for_this_experiment) && mkpath(result_folder_for_this_experiment)

    #the result file
    results_save_path = string(result_folder_for_this_experiment, "/ruin_depth_", now(), ".csv")

    df = DataFrame(depth=Float64[], Δ_fit=[], cpu_time=[])
    CSV.write(results_save_path, df)# write the header
 
    fill_adjacent_stations()
    #fill in the list of solutions
    rng = MersenneTwister(seed)
    solution_list = []
    for sc in scenario_list_ids
        initialize_scenarios([sc])
        push!(solution_list, [generate_random_solution() for _ in 1:trial_nbr])
    end

    for rd in ruin_depth_list
        @info "we are trying ruin depth = $rd"
        ruin_depth = rd
        sum_percentage = 0.0
        cpu_time = 0.0
        
        for sc in scenario_list_ids
            initialize_scenarios([sc])
            
            opt_path = "Data/MIP/solutions/E_carsharing_mip_scenario_$(sc)_requests_1000_walking_time_5.jls"
            fit = ECS_objective_function(load_sol(opt_path))
            
            for tr in 1:trial_nbr
                    sol = solution_list[sc][tr]
                    
                    sa_sol, sa_fit, sa_cpu = simulated_annealing(sol, 796.0, 1.0, 0.98, 82, 0.8)
                    
                    sum_percentage += (fit - sa_fit) / fit
                    cpu_time += sa_cpu
            end

        end
        empty!(df)
        push!(df, [ruin_depth, sum_percentage / (length(scenario_list) * trial_nbr), cpu_time / (length(scenario_list) * trial_nbr)])
        CSV.write(results_save_path, df, append=true)
    end
    
end

