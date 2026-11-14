using StatsBase
using Distributions
using DataFrames
using CSV
using ArgParse
include("model.jl")
include("evolution.jl")

user_args=ARGS[1:end]
name = user_args[1]
line1 = parse(Int,user_args[2])
line2 = parse(Int,user_args[3])

d = CSV.read("sims/"*name*".csv",DataFrame)
d = d[line1:line2,:]

if (length(user_args) > 3) && (user_args[4] == "Strategy3")
    type = Strategy3
else
    type = Strategy
end

pars = fieldnames(type)
for p in pars
    d[:,p] = [0.0 for _ in 1:nrow(d)]
end
d.fitness = [-1.0 for _ in 1:nrow(d)]

for i in 1:nrow(d)
    s = type(ntuple(_ -> 0.0, fieldcount(type))...)
    parameters = Dict("n" => d.n[i], "alpha" => d.alpha[i], "beta" => d.beta[i], "gamma" => d.gamma[i], "delta" => d.delta[i], "lambda" => d.lambda[i], "b" => d.b[i], "c" => d.c[i])
    fitnessFunction = x -> fitness(parameters, x)
    if d.mutationType[i] == "default"
        mutationFunction = x -> mutateStrategy(x, d.mutationProbability[i], d.mutationSize[i]; mutate_rho = false)
    elseif d.mutationType[i] == "all"
        mutationFunction = x -> mutateStrategy(x, d.mutationProbability[i], d.mutationSize[i]; mutate_rho = true)
    elseif d.mutationType[i] == "interactions_only"
        s = Strategy(pi = [1.0, 0.0, 0.0, 1.0], theta = [0.0, 0.0, 0.0, 0.0], rho = [0.0, 0.0])
        mutationFunction = x -> mutateStrategy(x, d.mutationProbability[i], d.mutationSize[i]; mutate_rho = false, mutate_pi = false)
    elseif d.mutationType[i] == "third_state_inactive"
        mutationFunction = x -> mutateStrategy(x, d.mutationProbability[i], d.mutationSize[i]; mutate_rho = true, third_state_active = false)
    end
    sim = evolutionarySimulation(s, fitnessFunction, strongSelectionFixationProbability, mutationFunction, d.nMutants[i])

    d.fitness[i] = sim[2]
    for p in pars
        d[i,p] = getproperty(sim[1],p)
    end
end

if (length(user_args) > 3) && (user_args[4] == "Strategy3")
    CSV.write("results/"*name*"_Strategy3/"*name*"_"*string(line1)*"_"*string(line2)*".csv",d)
else   
    CSV.write("results/"*name*"/"*name*"_"*string(line1)*"_"*string(line2)*".csv",d)
end

