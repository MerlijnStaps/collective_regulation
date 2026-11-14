push!(LOAD_PATH,"./")
push!(LOAD_PATH,"julia/")
include("model.jl")
using DataFrames, CSV

# A = good
# B = bad
strategies = Dict()
strategies["always_inactive"] = Strategy(0,1,0,1,0,0,0,0,0,0)
strategies["always_active"] = Strategy(1,0,1,0,0,0,0,0,0,0)
strategies["individual_regulation"] = Strategy(1,0,0,1,0,0,0,0,0,0)
strategies["collective_regulation_I"] = Strategy(1,0,0,1,1,0,0,0,0,0)
strategies["collective_regulation_II"] = Strategy(1,0,0,1,0,1,0,0,0,0)
strategies["collective_regulation_III"] = Strategy(1,0,0,1,0,0,0,1,0,0)
strategies["collective_regulation_IV"] = Strategy(1,0,0,1,0,0,1,0,0,0)

lambda = 10 .^ collect(-2:0.01:2)
beta = [0,1,4,10]
gamma = [0]
delta = [1e-6]
alpha = [1]
b = 2 .^collect(-4:0.02:4)
c = [1]
n = [0, 5]

d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,lambda,b,c,n)))
rename!(d,[:alpha,:beta,:gamma,:delta,:lambda,:b,:c,:n])
for (name,s) in strategies
    d[!, name * "_fitness"] = map(row ->fitness(Dict(String(k) => v for (k, v) in pairs(row)), s), eachrow(d))
end
d.winner = map(row -> collect(keys(strategies))[argmax([row[name * "_fitness"] for name in keys(strategies)])],eachrow(d))
d.winner_fitness = map(row -> maximum([row[name * "_fitness"] for name in keys(strategies)]),eachrow(d))
CSV.write("results/optimal_strategies.csv",d)