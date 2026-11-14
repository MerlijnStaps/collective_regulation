push!(LOAD_PATH,"./")
push!(LOAD_PATH,"julia/")
include("model.jl")
using Plots
using StatsPlots
using CSV
using DataFrames


for motif in ["I", "II", "spontaneous"]

    # Define strategy
    s = if motif == "I"
        Strategy(1, 0, 0, 1, 1, 0, 0, 0, 0, 0)
    elseif motif == "II"
        Strategy(1, 0, 0, 1, 0, 1, 0, 0, 0, 0)
    elseif motif == "spontaneous"
        Strategy(1, 0, 0, 1, 0, 0, 0, 0, sqrt(3), 0)
    end

    d = DataFrame(
        alpha = Float64[],
        beta = Float64[],
        x = Float64[],
        y = Float64[],
        env = String[]
    )

    for alpha in 1, beta in [0, 4]

        parameters = Dict(
            "alpha" => alpha,
            "beta" => beta,
            "lambda" => 1,
            "b" => 2,
            "c" => 1,
            "gamma" => 1,
            "delta" => 1e-6,
            "n" => 0
        )

        w = solveModel(parameters, s)
        x = w.I[1]:0.0001:w.I[2]

        append!(d, DataFrame(
            alpha = alpha,
            beta = beta,
            x = x,
            y = map(w.uA, x),
            env = "good"
        ))

        append!(d, DataFrame(
            alpha = alpha,
            beta = beta,
            x = x,
            y = map(w.uB, x),
            env = "bad"
        ))
    end

    CSV.write("results/motif_" * motif * "_examples.csv", d)

    # Generate mean values
    d_mean = DataFrame(
        alpha = Float64[],
        beta = Float64[],
        b = Float64[],
        meanA = Float64[],
        meanB = Float64[]
    )

    for alpha in 1, beta in [0,4], b in 2

        parameters = Dict(
            "alpha" => alpha,
            "beta" => beta,
            "lambda" => 1,
            "b" => b,
            "c" => 1,
            "gamma" => 1,
            "delta" => 1e-6,
            "n" => 0
        )

        w = solveModel(parameters, s)

        append!(d_mean, DataFrame(
            alpha = alpha,
            beta = beta,
            b = b,
            meanA = w.meanA,
            meanB = w.meanB
        ))
    end

    CSV.write("results/motif_" * motif * "_examples_means.csv", d_mean)
end
