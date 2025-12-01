#
#       Ambiente gráfico em Julia (muito melhor que o GnuPlot =V)
# ..................................................................................................
using CSV
using DataFrames
using GLMakie
using Printf
using ProgressMeter
# ..................................................................................................
#       Constantes.
#   → Pastas para ignorar na filtragem.
const folders_to_ignore = [".git", ".idea", "cmake-build-debug", "pr3c"]

#   → Arquivos para o painel de snapshots
const snapshots_id = [40, 80, 110, 130, 170, 240]
const snapshots_tag = [("(a)", [1, 1]), ("(b)", [1, 2]), ("(c)", [2, 1]), ("(d)", [2, 2]), ("(e)", [3, 1]), ("(f)", [3, 2])]

#   → Valores fixos
const const_gravitacional = 6.6743e-11
const mars_radius = 3.3895e6
const mars_mass = 6.4171e23
const m_to_ua = 1.496e11
# ..................................................................................................
#       Função principal.
function main()
    println("Executando script gráficos em '$(@__DIR__)' ...")
    # ..............................................................................................
    #       → Pega as pastas a fim de saber em qual delas executar o processamento.
    folders = filter(x -> isdir(joinpath(@__DIR__, x)), readdir(@__DIR__))
    folders = filter(x -> !(x in folders_to_ignore), folders)
    
    if (length(folders) < 1)
        println("Não foi identificado nenhum projeto. Por favor, rode a simulação primeiro.")
        return
    end
    
    println("Por favor, selecione um desses seguintes projetos para iniciar o processamento:")
    project::Int = 0

    for i in eachindex(folders)
        println("\t[$i] \t - Project: $(folders[i])")
    end

    while true
        print("Identificador do projeto: ")
        input = readline()
        project = parse(Int, input)
        if (project > 0 && project <= length(folders))
            break
        end
    end
    # ..............................................................................................
    print("Por favor, insira a velocidade no infinito (em m/s): ")
    velocidade_infinito = parse(Float64, readline())
    # ..............................................................................................
    mkpath(joinpath(@__DIR__, folders[project], "results"))
    # ..............................................................................................
    #       → Processa os dados globais.
    global_process(folders[project], velocidade_infinito)
    # ..............................................................................................
    #       → Gera snapshots da trajetória.
    snapshot_pr2c(folders[project])
    snapshot_pr3c(folders[project])
    # ..............................................................................................
end
# ..................................................................................................
#       Processamento dos dados por simulação.
function snapshot_pr2c(project_name::String)
    println("Inicializando processamento das trajetórias do projeto '$project_name' ...")
    # ..............................................................................................
    #       → Cria a pasta de processamento dos dados globais.
    input_path::String =  joinpath(@__DIR__, project_name, "pr2c")
    output_path::String = joinpath(@__DIR__, project_name, "results", "snapshots", "pr2c")
    mkpath(output_path)
    # ..............................................................................................
    #       → Lista os arquivos
    files = filter(x -> isfile(joinpath(input_path, x)), readdir(input_path))
    files = filter(x -> x != "global.csv", files)
    # ..............................................................................................
    #       → Pega os dados globais também, para a ter o parâmetro de impacto.
    data_global = CSV.read(joinpath(@__DIR__, project_name, "global_pr2c.csv"), DataFrame, delim = ",")
    data_global = Matrix(data_global)
    # ..............................................................................................
    #       → Painel de snapshots.
    painel = Figure(size = (1000, 1200))
    # ..............................................................................................
    #       → Abre cada um dos arquivos e processa eles
    @showprogress for file in files
        # ..........................................................................................
        #       → Lê o arquivo e pega os dados
        path = joinpath(input_path, file)
        data = CSV.read(path, DataFrame, delim = ",")
        data = Matrix(data)
        # ..........................................................................................
        #       → Pega o identificador
        m = match(r"\d+", file)
        i = parse(Int, m.match)
        # ..........................................................................................
        #       → Monta a figura
        fig = Figure(size = (600, 600))
        ax = Axis(fig[1, 1], title = @sprintf("Simulação para b = %.2e km", data_global[i, 2] / 1000), xlabel = L"x~[km]", ylabel = L"y~[km]")
        ax.xtickformat = "{:.2e}"
        ax.ytickformat = "{:.2e}"
        ax.xticks = Makie.LinearTicks(4)

        scatter!(ax, data[:, 2] ./ 1000, data[:, 3] ./ 1000, color = :blue, label = L"\text{Trajetória}", markersize = 4)

        scatter!(ax, 0.0, 0.0, markersize = 20, color = :red, label = L"\text{Marte}")
        axislegend(ax, position = :lc)
        # ..........................................................................................
        #       → Salva a figura.
        root, _ = splitext(file)
        save(joinpath(output_path, "$(root).png"), fig)
        # ..........................................................................................
        #       → Verifica para a montagem do painel.
        if i in snapshots_id
            tag, pos = snapshots_tag[findfirst(==(i), snapshots_id)]
            ax = Axis(painel[pos...], title = @sprintf("%s Simulação para b = %.2e km", tag, data_global[i, 2] / 1000), xlabel = L"x~[km]", ylabel = L"y~[km]")
            ax.xtickformat = "{:.2e}"
            ax.ytickformat = "{:.2e}"
            ax.xticks = Makie.LinearTicks(4)
            scatter!(ax, data[:, 2] ./ 1000, data[:, 3] ./ 1000, color = :blue, label = L"\text{Trajetória}", markersize = 4)
            scatter!(ax, 0.0, 0.0, markersize = 20, color = :red, label = L"\text{Marte}")
            axislegend(ax, position = :lc)
        end
    end
    # ..............................................................................................
    #       → Salva o painel.
    save(joinpath(@__DIR__, project_name, "results", "painel_pr2c.png"), painel)
    # ..............................................................................................
end
# ..................................................................................................
#       Processamento dos dados por simulação.
function snapshot_pr3c(project_name::String)
    println("Inicializando processamento das trajetórias do projeto '$project_name' ...")
    # ..............................................................................................
    #       → Cria a pasta de processamento dos dados globais.
    input_path::String =  joinpath(@__DIR__, project_name, "pr3c")
    output_path::String = joinpath(@__DIR__, project_name, "results", "snapshots", "pr3c")
    mkpath(output_path)
    # ..............................................................................................
    #       → Lista os arquivos
    files = filter(x -> isfile(joinpath(input_path, x)), readdir(input_path))
    files = filter(x -> x != "global.csv", files)
    # ..............................................................................................
    #       → Pega os dados globais também, para a ter o parâmetro de impacto.
    data_global = CSV.read(joinpath(@__DIR__, project_name, "global_pr3c.csv"), DataFrame, delim = ",")
    data_global = Matrix(data_global)
    # ..............................................................................................
    #       → Painel de snapshots.
    painel = Figure(size = (1000, 1200))
    # ..............................................................................................
    #       → Abre cada um dos arquivos e processa eles
    @showprogress for file in files
        # ..........................................................................................
        #       → Lê o arquivo e pega os dados
        path = joinpath(input_path, file)
        data = CSV.read(path, DataFrame, delim = ",")
        data = Matrix(data)
        # ..........................................................................................
        #       → Pega o identificador
        m = match(r"\d+", file)
        i = parse(Int, m.match)
        # ..........................................................................................
        #       → Monta a figura
        fig = Figure(size = (600, 600))
        ax = Axis(fig[1, 1], title = @sprintf("Simulação para b = %.2e km", data_global[i, 2] / 1000), xlabel = L"x~[U.A.]", ylabel = L"y~[U.A.]")
        ax.xtickformat = "{:.5f}"
        ax.ytickformat = "{:.2f}"
        ax.xticks = Makie.LinearTicks(3)

        scatter!(ax, data[:, 4] ./ m_to_ua, data[:, 5] ./ m_to_ua, color = :blue, label = L"\text{Trajetória da sonda}", markersize = 4)

        scatter!(ax, data[:, 2] ./ m_to_ua, data[:, 3] ./ m_to_ua, markersize = 4, color = :red, label = L"\text{Trajetória de Marte}")
        axislegend(ax, position = :cb)
        # ..........................................................................................
        #       → Salva a figura.
        root, _ = splitext(file)
        save(joinpath(output_path, "$(root).png"), fig)
        # ..........................................................................................
        #       → Verifica para a montagem do painel.
        if i in snapshots_id
            tag, pos = snapshots_tag[findfirst(==(i), snapshots_id)]
            ax = Axis(painel[pos...], title = @sprintf("%s Simulação para b = %.2e km", tag, data_global[i, 2] / 1000), xlabel = L"x~[U.A.]", ylabel = L"y~[U.A.]")
            ax.xtickformat = "{:.5f}"
            ax.ytickformat = "{:.2f}"
            ax.xticks = Makie.LinearTicks(3)
            scatter!(ax, data[:, 4] ./ m_to_ua, data[:, 5] ./ m_to_ua, color = :blue, label = L"\text{Trajetória da sonda}", markersize = 4)
            scatter!(ax, data[:, 2] ./ m_to_ua, data[:, 3] ./ m_to_ua, markersize = 4, color = :red, label = L"\text{Trajetória de Marte}")
            axislegend(ax, position = :cb)
        end
    end
    # ..............................................................................................
    #       → Salva o painel.
    save(joinpath(@__DIR__, project_name, "results", "painel_pr3c.png"), painel)
    # ..............................................................................................
end
# ..................................................................................................
#       Processamento dos dados globais.
function global_process(project_name::String, velocidade_infinito::Float64)
    println("Inicializando processamento dos dados globais do projeto '$project_name' ...")
    # ..............................................................................................
    #       → Cria a pasta de processamento dos dados globais.
    input_path::String =  joinpath(@__DIR__, project_name)
    output_path::String = joinpath(@__DIR__, project_name, "results", "global")
    mkpath(output_path)
    # ..............................................................................................
    #       → Lê o arquivo de dados globais.
    data_pr2c = CSV.read(joinpath(input_path, "global_pr2c.csv"), DataFrame, delim = ",")
    data_pr2c = Matrix(data_pr2c)

    data_pr3c = CSV.read(joinpath(input_path, "global_pr3c.csv"), DataFrame, delim = ",")
    data_pr3c = Matrix(data_pr3c)
    # ..............................................................................................
    #       Distância mínima por parâmetro de impacto.
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"\text{Parâmetro de impacto}~[km]", ylabel = L"\text{Distância mínima registrada}~[km]", title = "Distância mínima registrada pela sonda para cada parâmetro de impacto testado")

    d_esperado = zeros(Float64, length(data_pr2c[:, 2]))
    for i in eachindex(d_esperado)
        d_esperado[i] = - (const_gravitacional * mars_mass / velocidade_infinito^2) + sqrt((const_gravitacional * mars_mass / velocidade_infinito^2)^2 + data_pr2c[i, 2]^2)
    end

    colors_pr2c = ifelse.(data_pr2c[:, 6] .== 1, :red, :cyan)
    colors_pr3c = ifelse.(data_pr3c[:, 7] .== 1, :tomato, :black)

    scatter!(ax, data_pr2c[:, 2] ./ 1000, data_pr2c[:, 3] ./ 1000, markersize = 6, color = colors_pr2c)
    scatter!(ax, data_pr3c[:, 2] ./ 1000, data_pr3c[:, 3] ./ 1000, markersize = 8, color = (:blue, 0), marker = :rect, strokecolor = colors_pr3c, strokewidth = 1)
    hlines!(ax, mars_radius ./ 1000, label = L"R_{Marte}", color = :gray)

    scatter!(ax, [NaN], [NaN], color = :cyan, label = L"r_{min}~\text{(Trajetórias sem colisão, 2 corpos)}")
    scatter!(ax, [NaN], [NaN], color = :red, label = L"r_{min}~\text{(Trajetórias com colisão, 2 corpos)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :black, strokewidth = 1, label = L"r_{min}~\text{(Trajetórias sem colisão, 3 corpos)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :tomato, strokewidth = 1, label = L"r_{min}~\text{(Trajetórias com colisão, 3 corpos)}")

    lines!(ax, data_pr2c[:, 2] ./ 1000, d_esperado ./ 1000, color = :green, label = L"d_{min}~\text{(Valor teórico esperado, 2 corpos)}")

    axislegend(ax, position = :ct)

    save(joinpath(output_path, "distance.png"), fig)
    # ..............................................................................................
    #       Variação de velocidade por parâmetro de impacto.
    fig = Figure(size = (1720, 860))
    ax = Axis(fig[1, 1], xlabel = L"\text{Parâmetro de impacto}~[km]", ylabel = L"\text{Erro relativo}~[%]", title = "(a) Erro relativo da velocidade relativa por parâmetro de impacto simulado")

    scatter!(ax, data_pr2c[:, 2] ./ 1000, (data_pr2c[:, 4] ./ velocidade_infinito) .* 100, markersize = 6, color = colors_pr2c)

    v_init = sqrt(velocidade_infinito^2 + (2 * const_gravitacional * mars_mass) / (50 * mars_radius))
    data_nocol = []
    for i in eachindex(data_pr3c[:, 7])
        if data_pr3c[i, 7] == 1
            continue
        end
        
        push!(data_nocol, i)
    end

    scatter!(ax, data_pr3c[data_nocol, 2] ./ 1000, (abs.(data_pr3c[data_nocol, 5]) ./ v_init) .* 100, markersize = 8, color = (:blue, 0), marker = :rect, strokecolor = colors_pr3c[data_nocol], strokewidth = 1)

    scatter!(ax, [NaN], [NaN], color = :cyan, label = L"\Delta v/v_{∞}~\text{(Trajetórias sem colisão, 2 corpos)}")
    scatter!(ax, [NaN], [NaN], color = :red, label = L"\Delta v/v_{∞}~\text{(Trajetórias com colisão, 2 corpos)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :black, strokewidth = 1, label = L"\Delta v/|\vec{v}(0)|~\text{(Trajetórias sem colisão, 3 corpos)}")

    axislegend(ax, position = :lc)

    ax = Axis(fig[1, 2], xlabel = L"\text{Parâmetro de impacto}~[km]", ylabel = L"\text{Tempo de integração}~[h]", title = "(b) Tempo total de integração por parâmetro de impacto simulado")

    scatter!(ax, data_pr2c[:, 2] ./ 1000, data_pr2c[:, 7] ./ 3600, markersize = 6, color = colors_pr2c)
    scatter!(ax, data_pr3c[data_nocol, 2] ./ 1000, data_pr3c[data_nocol, 8] ./ 3600, markersize = 8, color = (:blue, 0), marker = :rect, strokecolor = colors_pr3c[data_nocol], strokewidth = 1)

    scatter!(ax, [NaN], [NaN], color = :blue, label = L"t~\text{(Trajetórias sem colisão)}")
    scatter!(ax, [NaN], [NaN], color = :red, label = L"t~\text{(Trajetórias com colisão)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :black, strokewidth = 1, label = L"t~\text{(Trajetórias sem colisão, 3 corpos)}")

    axislegend(ax, position = :rb)

    save(joinpath(output_path, "velocity.png"), fig)
    # ..............................................................................................
    #       Ângulo de deflexão.
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"\text{Parâmetro de impacto}~[km]", ylabel = L"\text{Ângulo de deflexão}~[∘]", title = "Ângulo de deflexão por parâmetro de impacto simulado")

    scatter!(ax, data_pr2c[:, 2] ./ 1000, data_pr2c[:, 5], markersize = 6, color = colors_pr2c)
    scatter!(ax, data_pr3c[:, 2] ./ 1000, data_pr3c[:, 6], markersize = 8, color = (:blue, 0), marker = :rect, strokecolor = colors_pr3c, strokewidth = 1)

    scatter!(ax, [NaN], [NaN], color = :blue, label = L"δ~\text{(Trajetórias sem colisão)}")
    scatter!(ax, [NaN], [NaN], color = :red, label = L"δ~\text{(Trajetórias com colisão)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :black, strokewidth = 1, label = L"δ~\text{(Trajetórias sem colisão, 3 corpos)}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :tomato, strokewidth = 1, label = L"δ~\text{(Trajetórias com colisão, 3 corpos)}")

    deflection_teo = zeros(Float64, length(data_pr2c[:, 2]))
    for i in eachindex(deflection_teo)
        deflection_teo[i] = rad2deg(2 * atan(const_gravitacional * mars_mass / (abs(data_pr2c[i, 2]) * velocidade_infinito^2)))
    end

    lines!(ax, data_pr2c[:, 2] ./ 1000, deflection_teo, color = :green, label = L"δ~\text{(Valor teórico esperado)}")

    axislegend(ax, position = :rt)

    save(joinpath(output_path, "deflection.png"), fig)
    # ..............................................................................................
    #       Variação da velocidade heliocêntrica da sonda
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"\text{Parâmetro de impacto}~[km]", ylabel = L"\text{Velocidade heliocêntrica}~[km/s]", title = "Variação da velocidade heliocêntrica da sonda")

    scatter!(ax, data_pr3c[:, 2] ./ 1000, data_pr3c[:, 4] ./ 1000, markersize = 8, color = (:blue, 0), marker = :rect, strokecolor = colors_pr3c, strokewidth = 1)

    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :black, strokewidth = 1, label = L"\Delta|\vec{v}|~\text{(Trajetórias sem colisão}")
    scatter!(ax, [NaN], [NaN], color = (:blue, 0), marker = :rect, strokecolor = :tomato, strokewidth = 1, label = L"\Delta|\vec{v}|~\text{(Trajetórias com colisão)}")
    axislegend(ax, position = :lt)

    save(joinpath(output_path, "velo_helio.png"), fig)
    # ..............................................................................................
end
# ..................................................................................................
#       Chamada da função principal.
out = main()