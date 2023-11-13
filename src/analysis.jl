"""
    plot_best_fitness_tracking(tracking_file_path::String; data_range::range=nothing)

    plot the progress of the best solution during the search process.
    @tracking_file_path: the path of the file containing the best fitness for each Iteration
    @data_range: the range of the data to be plotted (by default we plot all the data)
"""
function plot_best_fitness_tracking(tracking_file_path::String; data_range=nothing, optimal_fitness::Union{Nothing, Float64} = nothing)
    data_df = CSV.read(tracking_file_path, DataFrame, delim=';'; select=[:Iteration, :BestFitness], decimal=',')

    #delete all the points wher fit equal to 10 ^ 16
    #filter!(row -> row.BestFitness != 1e16, data_df)

    #select the range of data
    data_df = isnothing(data_range) ? data_df : data_df[data_range, :]
    #cleanup the plot 
    StatsPlots.plot();
    #plot the data
    if !isnothing(optimal_fitness)
        data_df.optimal_fitness .= optimal_fitness
        @df data_df StatsPlots.plot!(:Iteration, :optimal_fitness, linestyle=:dash , label = "optimal")
    end
    @df data_df StatsPlots.plot!(:Iteration, :BestFitness, xlabel="Iteration", ylabel="Fitness", title="Best Fitness")
    #ylims!(optimal_fitness - 10000, 0 )
end

"""
    read_heuristics_performance(performance_file_path::String; iteration_number::Int64)

    read the performance of the heuristics for a given iteration 
    it return a DataFrame of heuristicId the number of calls, the number of improvement,
    the number of equal and the number of worse solutions produced by heursitic.
    @performance_file_path: the path of the file containing the performance of the heuristics
    @iteration_number: the iteration number for which we want to read the performance 
                        (if isnothing it returns the last iteration)
    @save_file_path: the path of the file where the performance will be saved
"""
function read_heuristics_performance(performances_file_path::String; save_file_path::String = "", iter_num::Int64 = -1, data_range::Union{Nothing, UnitRange} = nothing)
    
    df = CSV.read(performances_file_path, DataFrame, delim = ";")
    
    # get the number of heuristics
    heuristics_number = Int64((size(df, 2) - 4) / 6)

    df = isnothing(data_range) ? df : df[data_range, :]

    # iteration number
    iter_num = (iter_num==-1) ? size(df, 1) : iter_num
    heuristics_performance_df = DataFrame(id = [], calls_number = [], best_found_number = [], improvement_number = [] , equal_number = [], worse_number = [], cpu_time = [])
    
    for i in 1:heuristics_number
        shift = 4+(i-1)*6
        push!(heuristics_performance_df,  (i, Tuple(df[iter_num, shift:(shift+5)])...) )
    end
    
    if save_file_path != ""
        CSV.write(save_file_path, heuristics_performance_df, decimal = ',')
    end
    sort!(heuristics_performance_df, :improvement_number, rev = true)
end

"""
    get_nbr_improvement(tracking_file_path::String)

see the the number of times that the solver could improve the overall best solution
"""
function get_nbr_improvement(tracking_file_path::String)
    data_df = CSV.read(tracking_file_path, DataFrame ,delim=';'; select = [:Iteration, :BestFitness], decimal = ',')
   
    counter = 0
    for i in 2:nrow(data_df)
        if data_df[i, :BestFitness] < data_df[i-1, :BestFitness]
            counter += 1
        end
    end
    @info "[solver Improvement capability]: $(counter) improvement(s) out of $(data_df[end, :Iteration]) iteration(s)"
    
    return counter
end

"""
    plot_preprocessing_results()
    plot the number of feasible trips for each  walking time number 
    of scenarios and number of requests per scenario
"""
function plot_preprocecessing_results()
    # Load data
    df = CSV.read("results/preprocessing_2019/preprocessing_results.csv", DataFrame)
    # the data that we will group on to reshape our data
    columns_group = [:K, :β_w]
    data_as_groups = groupby(df, columns_group)
    #make data as Matrix to plot the groupedbars
    data_as_Matrix = Array{Int,2}[]
    for gr in data_as_groups
        if isempty(data_as_Matrix)
            data_as_Matrix = gr.H'
        else
            data_as_Matrix = vcat(data_as_Matrix, gr.H')
        end
    end


    keys(data_as_groups)
    layout = @layout [a b; c d e; f{0.001h}]
    titles = ["(a): scenarios of 1000 requests",
        "(b): scenarios of 2000 requests",
        "(c): scenarios of 3000 requests",
        "(d): scenarios of 4000 requests",
        "(e): scenarios of 5000 requests",
        "The total number of feasible trips genrated"]

    p_objs = []
    for i in 1:5

        data = data_as_Matrix[((i-1)*7+1):(i*7), :]
        p = groupedbar(data,
            title="$(titles[i]): ",
            bar_width=1.0,
            bar_position=:dodge, legend=false,
            ylabel = "Number of feasible trips"
        )
        push!(p_objs, p)
    end
    #add the title plots
    push!(p_objs, plot(legend=false, showaxis=false, grid=false, title=titles[6]))

    p = plot(p_objs..., layout=layout, size=(1800, 1200))

end

"""
    plot_100_500_scenario_results()

"""
function plot_100_500()
    #= ***************** the first plot of feasile trips for 1000 requests per scenario =#
    # Load data
    df = CSV.read("results/preprocessing_2019/preprocessing_results.csv", DataFrame)
    # the data that we will group on to reshape our data
    columns_group = [:K, :β_w]
    data_as_groups = groupby(df, columns_group)
    #make data as Matrix to plot the groupedbars
    data_as_Matrix = Array{Int,2}[]
    for gr in data_as_groups
        if isempty(data_as_Matrix)
            data_as_Matrix = gr.H'
        else
            data_as_Matrix = vcat(data_as_Matrix, gr.H')
        end
    end

    data_1000_requests_per_sceanrio = data_as_Matrix[1:7, :]
    p1 = groupedbar(data_1000_requests_per_sceanrio,
        bar_width=1.0,
        bar_position=:dodge, label=["100 scenario" "200 scenario" "300 scenario" "400 scenario" "500 scenario"],
        ylabel="|H|",
        xlabel="β_w")

    x_custom_tickets = ["5", "6", "7", "8", "9", "10", "15"]
    xticks!(p1, (1:7, x_custom_tickets))
    title!("Feasible trips generated for 1000 requests per scenario", titlefont=font(10))
    #= ********** the second plot for the feasible paths *************** =#
    data_5000_requests_per_sceanrio = data_as_Matrix[29:35, :]
    p2 = groupedbar(data_5000_requests_per_sceanrio,
        bar_width=1.0,
        bar_position=:dodge, label=["100 scenario" "200 scenario" "300 scenario" "400 scenario" "500 scenario"],
        ylabel="|H|",
        xlabel="β_w")

    x_custom_tickets = ["5", "6", "7", "8", "9", "10", "15"]
    xticks!(p2, (1:7, x_custom_tickets))
    title!("Feasible trips generated for 5000 requests per scenario", titlefont=font(10))
    #= ************* plots the number of feasibles requests ************************** =#
    data_as_Matrix = Array{Int,2}[]
    for gr in data_as_groups
        if isempty(data_as_Matrix)
            data_as_Matrix = gr.K_a'
        else
            data_as_Matrix = vcat(data_as_Matrix, gr.K_a')
        end
    end

    data_1000_requests_per_sceanrio = data_as_Matrix[1:7, :]
    p3 = groupedbar(data_1000_requests_per_sceanrio,
        bar_width=1.0,
        bar_position=:dodge, label=["100 scenario" "200 scenario" "300 scenario" "400 scenario" "500 scenario"],
        ylabel="|Kₐ|",
        xlabel="β_w")

    x_custom_tickets = ["5", "6", "7", "8", "9", "10", "15"]
    xticks!(p3, (1:7, x_custom_tickets))
    title!("Feasible requests for 1000 requests per scenario", titlefont=font(10))
    #= ********** the second plot for the feasible paths *************** =#
    data_5000_requests_per_sceanrio = data_as_Matrix[29:35, :]
    p4 = groupedbar(data_5000_requests_per_sceanrio,
        bar_width=1.0,
        bar_position=:dodge, label=["100 scenario" "200 scenario" "300 scenario" "400 scenario" "500 scenario"],
        ylabel="|Kₐ|",
        xlabel="β_w")

    x_custom_tickets = ["5", "6", "7", "8", "9", "10", "15"]
    xticks!(p4, (1:7, x_custom_tickets))
    title!("Feasible requests for 5000 requests per scenario", titlefont=font(10))
    #save the different plots
    savefig(p1, "results/preprocessing_2019/feasible_trips_1000_requests_per_scenario.svg")
    savefig(p2, "results/preprocessing_2019/feasible_trips_5000_requests_per_scenario.svg")
    savefig(p3, "results/preprocessing_2019/feasible_requests_1000_requests_per_scenario.svg")
    savefig(p4, "results/preprocessing_2019/feasible_requests_5000_requests_per_scenario.svg")
end

function plot_heuristics_performance(performances_file_path::String;
                                             iter_num::Int64 = -1,
                                             data_range::Union{Nothing, UnitRange} = nothing)
    df = CSV.read(performances_file_path, DataFrame, delim = ";", decimal=',')
    
    # get the number of heuristics
    heuristics_number = Int64((size(df, 2) - 4) / 6)
    
    #get the data for the given range if it is set
    df = isnothing(data_range) ? df : filter(x-> x.Iteration ∈ data_range, df)
    
    heuristics_performance_df = DataFrame(id = [], calls_number = [], best_found_number = [], improvement_number = [] , equal_number = [], worse_number = [], cpu_time = [])
    
    for i in 1:heuristics_number
        shift = 4+(i-1)*6
        curr_heur_data = Vector(df[end, shift:(shift+5)]) - Vector(df[1, shift:(shift+5)]) 
        push!(heuristics_performance_df,  (i, curr_heur_data...) )
    end
    
    #sort!(heuristics_performance_df, focus_on, rev = true)

    positions = 1:size(heuristics_performance_df, 1)
    ticks = "LHH-" .* string.(heuristics_performance_df[:, :id] .- 1)
    @df heuristics_performance_df bar(:improvement_number, label="Improvement solutions found" , legendposition = :topright)
    @df heuristics_performance_df bar!(:best_found_number, label="new best solutions found")
    StatsPlots.xticks!((collect(positions), ticks))
    
end

"""
    get_best_fitness()
    return the fitness value of the best solution found by the solver

"""
function get_best_fitness(tracking_file_path::String)
    data_df = CSV.read(tracking_file_path, DataFrame ,delim=';'; #= select = [:Iteration, :BestFitness], =# decimal = ',')
   
    return data_df[end, :BestFitness]
end


"""
    plot_depth_cpu_fitness_experiment(file_path::String)

    plot two scatter plots of the fitness and cpu time for each initial_nmber of initial depth of search.
    the two scatters plots are plotted on the same figure after normalizing the values.
    The resulted plots allows as to select the best initial depth value
"""
function plot_depth_cpu_fitness_experiment(file_path::String)
    
    df = CSV.read(file_path, DataFrame)
    df.mean_fitness .*= -1 
    a = -1 .* df.cpu_time .+ maximum(df.cpu_time) .+ 1

    #normalize the fitness andc cpu
    df.y_values = df.mean_fitness ./ maximum(df.mean_fitness) 
    df.z_values = a ./ maximum(a)
    df.x_values = string.(df.init_depth, "_", df.iterations)

    #plot line plot
    scatter(df.x_values, df.y_values, label="Fit")
    scatter!(df.x_values, df.z_values, xrotation = 90, label= "CPU")
    xticks!(1:length(df.x_values), df.x_values)

end