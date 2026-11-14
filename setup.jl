using StatsBase
using Distributions
using DataFrames
using CSV
using Random

alpha = 10 .^ collect(-1:0.025:2)
beta = [1.0, 4.0, 10.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
gamma = [10.0]
delta = [1e-6]
lambda = [1.0]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/N1.csv",d)



alpha = 10 .^ collect(-1:0.1:2)
beta = [4.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
gamma = [10.0]
delta = [1e-6]
lambda = [1.0]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all","third_state_inactive"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/M1old.csv",d)




lambda = 10 .^ collect(-1:0.02:1)
beta = [1.0, 4.0, 10.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
alpha = [1.0]
gamma = [10.0]
delta = [1e-6]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/L1.csv",d)





lambda = 10 .^ collect(-1:0.04:1)
beta = [40.0, 100.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
alpha = [1.0]
gamma = [10.0]
delta = [1e-6]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/L2.csv",d)


lambda = 10 .^ collect(-1:0.04:1)
beta = [4.0, 40.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
alpha = [1.0]
gamma = [100.0]
delta = [1e-6]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/L3.csv",d)




lambda = 10 .^ collect(-1:0.05:1)
beta = [4.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
gamma = [10.0]
delta = [1e-6]
alpha = [1.0]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all","third_state_inactive"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/M1.csv",d)



lambda = 10 .^ collect(-1:0.04:1)
beta = [4.0,40.0]
b = 2 .^ collect(-3:0.1:3)
c = [1.0]
gamma = [100.0]
delta = [1e-6]
alpha = [1.0]
n = [5]
replicate = [1,2,3,4,5,6,7,8,9,10]
mutationType = ["all","third_state_inactive"]
d = DataFrame(collect(Iterators.product(alpha,beta,gamma,delta,b,c,lambda,n,replicate,mutationType)))
rename!(d,[:alpha,:beta,:gamma,:delta,:b,:c,:lambda,:n,:replicate,:mutationType])
d.nMutants = [100000 for _ in 1:nrow(d)]
d.mutationSize = [0.3 for _ in 1:nrow(d)]
d.mutationProbability = [0.1 for _ in 1:nrow(d)]
Random.seed!(123)
d = d[shuffle(1:nrow(d)), :]
CSV.write("sims/M2.csv",d)