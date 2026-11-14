push!(LOAD_PATH,"./")
push!(LOAD_PATH,"julia/")
include("model.jl")
using Plots
using DataFrames
using CSV

function parameters(row)
    Dict("alpha" => row["alpha"], "beta" => row["beta"], "gamma" => row["gamma"], "delta" => row["delta"], "lambda" => row["lambda"], "b" => row["b"], "c" => row["c"], "n" => row["n"])
end
function strategy(row)
    Strategy(row["piAA"],row["piAB"],row["piBA"],row["piBB"],row["thetaAA"],row["thetaAB"],row["thetaBA"],row["thetaBB"],row["rhoA"],row["rhoB"])
end
function process_results(name)
    d1 = CSV.read("results/"*name*"_results.csv.gz",DataFrame)
    d2 = d1[(d1.n .== 5), :]
    results = [solveModel(parameters(row),strategy(row)) for row in eachrow(d2)]

    d2.meanA = [clamp(x.meanA,0,1) for x in results]
    d2.meanB = [clamp(x.meanB,0,1) for x in results]
    for i in 1:6
        d2[!, "uA_" * string(i-1)] = [clamp(x.uA[i],0,1) for x in results]
        d2[!, "uB_" * string(i-1)] = [clamp(x.uB[i],0,1) for x in results]
    end
    CSV.write("results/"*name*"_processed_new.csv", d2)
end
process_results("L3")