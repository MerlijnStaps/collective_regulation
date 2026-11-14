push!(LOAD_PATH,"./")
push!(LOAD_PATH,"julia/")
include("model.jl")
using DataFrames, CSV

gamma = [0]
delta = [1e-6]
b = [2]
c = [1]
n = [0]

strategies = Dict()
strategies["always_inactive"] = Strategy(0,1,0,1,0,0,0,0,0,0)
strategies["always_active"] = Strategy(1,0,1,0,0,0,0,0,0,0)
strategies["individual_regulation"] = Strategy(1,0,0,1,0,0,0,0,0,0)
strategies["collective_regulation_I"] = Strategy(1,0,0,1,1,0,0,0,0,0)
strategies["collective_regulation_II"] = Strategy(1,0,0,1,0,1,0,0,0,0)
strategies["collective_regulation_III"] = Strategy(1,0,0,1,0,0,0,1,0,0)
strategies["collective_regulation_IV"] = Strategy(1,0,0,1,0,0,1,0,0,0)

alpha = 10 .^ collect(-2:0.02:2)
beta = collect(0:0.1:12)
lambda = [1]

d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,lambda,b,c,n)))
rename!(d,[:alpha,:beta,:gamma,:delta,:lambda,:b,:c,:n])
for (name,s) in strategies
    d[!, name * "_fitness"] = map(row ->fitness(Dict(String(k) => v for (k, v) in pairs(row)), s), eachrow(d))
end
d.winner = map(row -> collect(keys(strategies))[argmax([row[name * "_fitness"] for name in keys(strategies)])],eachrow(d))
d.winner_fitness = map(row -> maximum([row[name * "_fitness"] for name in keys(strategies)]),eachrow(d))

CSV.write("results/optimal_interaction_rates_alpha.csv",d) # takes 1 mins

alpha = [1,3,5]
beta = collect(0:0.05:12)
lambda = 10 .^ collect(-2:0.05:2)

d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,lambda,b,c,n)))
rename!(d,[:alpha,:beta,:gamma,:delta,:lambda,:b,:c,:n])
for (name,s) in strategies
    d[!, name * "_fitness"] = map(row ->fitness(Dict(String(k) => v for (k, v) in pairs(row)), s), eachrow(d))
end
d.winner = map(row -> collect(keys(strategies))[argmax([row[name * "_fitness"] for name in keys(strategies)])],eachrow(d))
d.winner_fitness = map(row -> maximum([row[name * "_fitness"] for name in keys(strategies)]),eachrow(d))

CSV.write("results/optimal_interaction_rates_lambda.csv",d) # takes 2 mins

alpha = collect(0.002:0.002:1)
lambda = 10 .^ collect(-2:0.05:2)
k = [4]
alpha0 = [1]

d = DataFrame(collect(Iterators.product(alpha,gamma,delta,lambda,b,c,n,k,alpha0)))
rename!(d,[:alpha,:gamma,:delta,:lambda,:b,:c,:n,:k,:alpha0])
d.beta = d.k .* (d.alpha0 .- d.alpha)
for (name,s) in strategies
    d[!, name * "_fitness"] = map(row -> fitness(Dict(String(k) => v for (k, v) in pairs(row)), s), eachrow(d))
end
d.winner = map(row -> collect(keys(strategies))[argmax([row[name * "_fitness"] for name in keys(strategies)])],eachrow(d))
d.winner_fitness = map(row -> maximum([row[name * "_fitness"] for name in keys(strategies)]),eachrow(d))

CSV.write("results/optimal_interaction_rates_costs.csv",d) # takes 1 mins