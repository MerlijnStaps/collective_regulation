
function mutateStrategy(s, mutationProbability, mutationSize; mutate_pi = true, mutate_theta = true, mutate_rho = true, third_state_active = true)
    pars = []
    d = Dict()
    for name in fieldnames(typeof(s))
        value = getfield(s, name)
        d[name] = value
    end
    if mutate_pi
        append!(pars, filter(name -> startswith(String(name), "pi"), fieldnames(typeof(s))))
    end
    if mutate_theta
        append!(pars, filter(name -> startswith(String(name), "theta"), fieldnames(typeof(s))))
    end
    if mutate_rho
        append!(pars, filter(name -> startswith(String(name), "rho"), fieldnames(typeof(s))))
    end
    for par in pars
        mutation_sign = rand() < 1/3 ? 1 : -1
        mutation = mutationSize * mutation_sign * abs(randn())
        if rand() < mutationProbability
            d[par] = clamp(d[par] + mutation, 0.0, 1.0) 
        end
    end
    for par in filter(name -> startswith(String(name), "activity"), fieldnames(typeof(s)))
        if third_state_active
            d[par] = 1.0
        else
            d[par] = 0.0
        end
    end
    # reduce thetas
    if typeof(s) == Strategy
        d[:thetaAA],d[:thetaBA] = reduce_pair(d[:thetaAA],d[:thetaBA])
        d[:thetaAB],d[:thetaBB] = reduce_pair(d[:thetaAB],d[:thetaBB])
    end
    if typeof(s) == Strategy3
        d[:thetaXXY], d[:thetaYXY], d[:thetaZXY] = reduce_triple(d[:thetaXXY], d[:thetaYXY], d[:thetaZXY])
        d[:thetaXXZ], d[:thetaYXZ], d[:thetaZXZ] = reduce_triple(d[:thetaXXZ], d[:thetaYXZ], d[:thetaZXZ])
        d[:thetaXYX], d[:thetaYYX], d[:thetaZYX] = reduce_triple(d[:thetaXYX], d[:thetaYYX], d[:thetaZYX])
        d[:thetaXYZ], d[:thetaYYZ], d[:thetaZYZ] = reduce_triple(d[:thetaXYZ], d[:thetaYYZ], d[:thetaZYZ])
        d[:thetaXZX], d[:thetaYZX], d[:thetaZZX] = reduce_triple(d[:thetaXZX], d[:thetaYZX], d[:thetaZZX])
        d[:thetaXZY], d[:thetaYZY], d[:thetaZZY] = reduce_triple(d[:thetaXZY], d[:thetaYZY], d[:thetaZZY])
    end
    # normalize pis
    for prefix in ["piAX","piAY","piAZ","piBX","piBY","piBZ","thetaXX","thetaXY","thetaXZ","thetaYX","thetaYY","thetaYZ","thetaZX","thetaZY","thetaZZ"]
        pars = filter(name -> startswith(String(name), prefix), fieldnames(typeof(s)))
        if length(pars) > 0
            total = sum(d[k] for k in pars)
            if total > 1
                for k in pars
                    d[k] /= total
                end
            end
        end
    end
    T = typeof(s)
    return T([d[Symbol(f)] for f in fieldnames(T)]...)
end

function mutateValue(x, mutationProbability, mutationSize, p_positive = 0.5)
    if rand() < mutationProbability
        mutation_sign = rand() < p_positive ? 1 : -1
        x = clamp(x + mutationSize * mutation_sign * abs(randn()), 0.0, 1.0)
    end
    return x
end

function reduce_pair(a,b)
    if a > 0 && b > 0
        m = min(a,b)
        a -= m
        b -= m
    end
    return a, b
end


function reduce_triple(a,b,c)
    if a > 0 && b > 0 && c > 0
        m = min(a,b,c)
        a -= m
        b -= m
        c -= m
    end
    return a, b, c
end

function evolutionarySimulation(x, fitness, fixation, mutation, n; verbose=true)
    f = fitness(x)
    if verbose
        println("The initial fitness is " * string(f))
    end
    count = 0
    failcount = 0
    for _ in 1:n
        y = mutation(deepcopy(x))
        #println(y)
        g = fitness(y)
        fixationProbability = fixation(g, f) 
        if rand() < fixationProbability
            if verbose 
                println("After " * string(failcount) * " mutants failed, " * string(x) * " fixed with fitness = " * string(g) * " and fixation probability " * string(fixationProbability))
            end
            x = y
            f = g
            count += 1
            failcount = 0
        else
            failcount += 1
        end
    end
    if verbose
        println(string(failcount) *" more mutants failed")
    end
    println(string(count) * " mutants fixed")
    return x, f
end

function strongSelectionFixationProbability(mutantFitness, residentFitness)
    return mutantFitness > residentFitness ? 1 : 0
end
