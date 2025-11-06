using CSV
using GLMakie
using DataFrames
using ProgressMeter

r_m = 2.28e11   # raio fixo da órbita de Marte
ua  = 1.496e11  # valor de 1 unidade atômica.

global_data = Matrix(CSV.read("data/global.dat", DataFrame, delim = " ", header = false))

b           = global_data[:, 1] ./ ua
delta_v     = global_data[:, 2] ./ 1000
deflexao    = global_data[:, 3]

g_fig = Figure(size = (1600, 800))
g_ax_delta = Axis(g_fig[1, 1], xlabel = "b [U.A.]", ylabel = "Δv [km/s]")
g_ax_def = Axis(g_fig[1, 2], xlabel = "b [U.A.]", ylabel = "δ [rad]")
scatter!(g_ax_delta, b, delta_v, markersize = 8)
scatter!(g_ax_def, b, deflexao, markersize = 8)

save("global.png", g_fig)



data_polar = Matrix(CSV.read("data/result_30.dat", DataFrame, delim = " ", header = false))


#       Coordenadas polares.
t   = data_polar[:, 1]
θ_m = data_polar[:, 2]
r_s = data_polar[:, 3]
θ_s = data_polar[:, 4]



#       Coordenadas cartesianas.
x_m = r_m .* cos.(θ_m)
y_m = r_m .* sin.(θ_m)

x_s = r_s .* cos.(θ_s)
y_s = r_s .* sin.(θ_s)

x_m ./= ua
y_m ./= ua
x_s ./= ua
y_s ./= ua

#       Build the figure
time = Observable(t[1])
points_ship = Observable(Point2f[(x_s[1], y_s[1])])

x_ship = Observable(x_s[1])
y_ship = Observable(y_s[1])

x_mars = Observable(x_m[1])
y_mars = Observable(y_m[1])

lims_obs = lift(x_mars, y_mars) do xm, ym
    (xm - 0.015, xm + 0.015, ym - 0.015, ym + 0.015)
end

fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1], xlabel = L"x~[U.A.]", ylabel = L"y~[U.A.]", title = lift(f -> "Time: $(round(f, digits=2)) h", time))
scatter!(ax, points_ship, markersize = 2, marker = :rect, color = :gray)
scatter!(ax, x_mars, y_mars, marker = :circle, markersize = 20, color = :red)
scatter!(ax, x_ship, y_ship, marker = '✈', markersize = 10, color = :blue)

on(lims_obs) do (xmin, xmax, ymin, ymax)
    limits!(ax, xmin, xmax, ymin, ymax)
end

frames = 2:length(t)

record(fig, "ship.mp4", frames; framerate = 30) do f
    new_point = Point2f(x_s[f], y_s[f])
    points_ship[] = push!(points_ship[], new_point)
    x_mars[] = x_m[f]
    y_mars[] = y_m[f]
    x_ship[] = x_s[f]
    y_ship[] = y_s[f]

    time[] = t[f] / 3600
end
