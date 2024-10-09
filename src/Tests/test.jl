using CSV, DataFrames
using Plots, StatsPlots
using Colors
using Statistics
mm = Plots.mm

df = CSV.read("results/solve_single_scenario_with_SA/2024-03-15T09:28:17.511/SA_ruin_recreate_details2024-03-15T09:28:17.511.csv", DataFrame)

data = Matrix{Float64}(undef, 10, 10)
for i in 1:10
    curr_sc_data = df[df.scenario_id .== i, :gap]
    data[:, i] = curr_sc_data
end

red_color = RGB(134/255, 150/255, 166/255)   # Custom red

# Define alternating colors
colors = [red_color, red_color]  # Alternating blue and red
bockground = RGB(225/255, 225/255, 225/255)  # Custom background color
# Initialize the plot
plt = plot(legend = false, xlabel = "Scenario ID", ylabel = "Gap", 
           title = "Gap distribution for each scenario", background_color = bockground,
           size = (1400, 600),
        left_margin = 10mm,
        bottom_margin= 10mm) # Adjust left and right margin

# Add each series to the boxplot with alternating custom colors
for i in 1:size(data, 2)
    # Shift the position of the box plot slightly using `i` to prevent overlapping
    boxplot!(i * ones(size(data, 1)), data[:, i], 
             color = colors[i % 2 + 1], 
             linecolor = :black, 
             markerstrokecolor = :black)
end

# Display the plot
display(plt)

all_data = data[:]  
boxplot!(11 * ones(size(all_data, 1)), all_data, color = colors[1],linecolor = :black, 
             markerstrokecolor = :black)

x_ticks = vcat(["sc$i" for i in 1:10], ["All"])
xticks!(1:11, x_ticks) 
y_ticks = ["0%", "0.5%", "1%", "1.5%"]
yticks!([0, 0.5, 1, 1.5], y_ticks)


agg_df = CSV.read("results/solve_single_scenario_with_SA/2024-03-15T09:28:17.511/SA_ruin_recreate_aggregated2024-03-15T09:28:17.511.csv", DataFrame)

y = -1 .* agg_df[1:10, :mean_gap]

rand_agg_df = CSV.read("results/solve_single_scenario_with_SA/2024-03-20T09:51:37.665/SA_ruin_recreate_aggregated2024-03-20T09:51:37.665.csv", DataFrame)
y2 = -1 .*rand_agg_df[1:10, :mean_gap]


plot(1:10, y, label = "SA + Adjacent Ruin", xlabel = "Scenario ID", ylabel = "Gap", 
title = "Average Gaps", background_color = bockground,
size = (1400, 600),
left_margin = 10mm,
bottom_margin= 10mm, legend=:right)
scatter!(1:10, y,label="")
plot!(1:10, y2, label = "SA + Random Ruin")
scatter!(1:10, y2, label="")

xticks!(1:10, x_ticks[1:10])
yticks!([0, 1, 2, 4, 6, 8], ["0%", "1%", "2%", "4%", "6%", "8%"])

scatter([0.0, 0.3],[0.0, 2.8], background_color=:transparent)

xticks!(collect(0.00:0.01:0.3), string.(collect(0.00:0.01:0.3)))



yticks!(collect(0.0:0.1:2.5), string.(collect(0.0:0.1:2.5)))
savefig("test.png")

    y = [2.45, 0.95, 0.85 , 0.97, 0.75, 0.73, 1.06, 1.08, 1.18, 0.95, 1.05, 0.93, 1.38, 1.34, 1.20, 1.34, 1.19, 1.15, 1.22, 1.36, 1.35, 1.120, 1.32, 1.25, 1.05, 1.45, 1.35, 1.39, 1.8, 1.98, 1.97]

    x = collect(0:0.01:(length(y)-1)*0.01)
    plot(x, y, size=(1400, 600), background_color = bockground, legend=false, xlabel="γ", ylabel="Gap",left_margin = 10mm,
    bottom_margin= 10mm, title="Aggregated Gaps for different γ values")
    scatter!(x, y)
    xticks!([0., 0.1, 0.2, 0.3], string.([0., 0.1, 0.2, 0.3]))
    yticks!([1.0, 1.5, 2.0, 2.5], ["$i%" for i in [ 1.0, 1.5, 2.0, 2.5]])

