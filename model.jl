using QuadGK
using LinearAlgebra

# A is good, B is bad

struct Strategy
    piAA::Float64 # probability to become active in good conditions 
    piAB::Float64
    piBA::Float64
    piBB::Float64
    thetaAA::Float64 # probability to become active upon meeting an active individual
    thetaAB::Float64
    thetaBA::Float64
    thetaBB::Float64
    rhoA::Float64 # probability to become active spontaneously
    rhoB::Float64 
end

Strategy(; pi::Vector{Float64},
         theta::Vector{Float64},
         rho::Vector{Float64}) = begin

    length(pi) == 4 || error("pi must have length 3")
    length(theta) == 4 || error("theta must have length 3")
    length(rho) == 2 || error("rho must have length 3")

    Strategy(
        pi[1], pi[2], pi[3], pi[4],
        theta[1], theta[2], theta[3], theta[4],
        rho[1], rho[2]
    )
end

function solveModel(parameters, s)
    if typeof(s) == Strategy
        if parameters["n"] == 0
            sol = solveInfiniteModel(parameters, s)
            I = sol.I
        else
            corrected_parameters = copy(parameters)
            corrected_parameters["beta"] = corrected_parameters["beta"] / (corrected_parameters["n"] - 1)
            sol = solveFiniteModel(corrected_parameters, s)
            I = [0.0, 1.0]
        end
    elseif typeof(s) == Strategy3
            corrected_parameters = copy(parameters)
            corrected_parameters["beta"] = corrected_parameters["beta"] / (corrected_parameters["n"] - 1)
            sol = solveFiniteModel3(corrected_parameters, s)
        I = [0.0, 1.0]
    end
    return (uA = sol.uA, uB = sol.uB, meanA = sol.meanA, meanB = sol.meanB, I = I) 
end

function fitness(parameters, s)
    sol = solveModel(parameters, s)
    return parameters["b"]*sol.meanA - parameters["c"]*sol.meanB
end

function solveFiniteModel(parameters,s)
    m = createSystemMarkovChain(parameters["n"], parameters["alpha"], parameters["beta"], parameters["gamma"], parameters["delta"], parameters["lambda"], s)
    n = parameters["n"]
    v = stationaryState(m)    
    uA = v[1:(n+1)]
    uB = v[(n+2):(2*n+2)] 
    uA = uA/sum(uA)
    uB = uB/sum(uB)
    meanA = sum(uA .* (0:n))/n
    meanB = sum(uB .* (0:n))/n 
    return(uA=uA,uB=uB,meanA=meanA,meanB=meanB)
end

function createSystemMarkovChain(n, alpha, beta, gamma, delta, lambda, s)
    m = zeros(2*n+2, 2*n+2)
    for i in 1:(n+1)
        m[i,i] += -lambda
        m[n+1+i, i] += lambda
        m[n+1+i, n+1+i] += -lambda
        m[i, n+1+i] += lambda
    end
    function pIncrease(x,env) # probability for A to increase given x individuals in env
        if env == "A"
            return alpha*(n-x)*s.piAA + beta*x*(n-x)*s.thetaAA + beta*(n-x)*(n-x-1)*s.thetaBA + (n-x)*gamma*s.rhoA + (n-x)*delta 
        elseif env == "B"
            return alpha*(n-x)*s.piBA + beta*x*(n-x)*s.thetaAA + beta*(n-x)*(n-x-1)*s.thetaBA + (n-x)*gamma*s.rhoA + (n-x)*delta
        end
    end
    function pDecrease(x, env)
        if env == "A"
            return alpha*x*s.piAB + beta*x*(n-x)*s.thetaBB + beta*x*(x-1)*s.thetaAB + x*gamma*s.rhoB + x*delta
        elseif env == "B"
            return alpha*x*s.piBB + beta*x*(n-x)*s.thetaBB + beta*x*(x-1)*s.thetaAB + x*gamma*s.rhoB + x*delta
        end    
    end
    for i in 1:(n)
        m[i,i+1] += pDecrease(i, "A")
        m[i+1, i+1] += -pDecrease(i, "A")
        m[i+1,i] += pIncrease(i-1, "A")
        m[i,i] += -pIncrease(i-1, "A")
        m[n+1+i,n+1+i+1] += pDecrease(i, "B")
        m[n+1+i+1,n+1+i+1] += -pDecrease(i, "B")
        m[n+1+i+1,n+1+i] += pIncrease(i-1, "B")
        m[n+1+i,n+1+i] += -pIncrease(i-1, "B")
    end
    return m
end


function solveInfiniteModel(parameters, s)
    alpha, beta, gamma, delta, lambda = parameters["alpha"], parameters["beta"], parameters["gamma"], parameters["delta"], parameters["lambda"] 

    function hA(x)
        return(alpha*s.piAA*(1-x) - alpha*s.piAB*x + beta*x*(1-x)*(s.thetaAA-s.thetaBB) + beta*(1-x)*(1-x)*s.thetaBA - beta*x*x*s.thetaAB + gamma*(1-x)*s.rhoA - gamma*x*s.rhoB + delta*(0.5-x))
    end
    function hB(x)
        return(alpha*s.piBA*(1-x) - alpha*s.piBB*x + beta*x*(1-x)*(s.thetaAA-s.thetaBB) + beta*(1-x)*(1-x)*s.thetaBA - beta*x*x*s.thetaAB + gamma*(1-x)*s.rhoA - gamma*x*s.rhoB + delta*(0.5-x))
    end

    # if hA = 0 or hB = 0, then the dynamics are stationary in one environment
    if hA(0) == hA(0.5) == hA(1) == 0 && hB(0) == hB(0.5) == hB(1) == 0
        return(hA = hA, hB = hB, meanA = 0.5, meanB = 0.5, uA = x -> 0.0, uB = x -> 0.0)
    elseif hA(0) == hA(0.5) == hA(1) == 0
        return(hA = hA, hB = hB, meanA = quadraticRoots(hB)[1].value, meanB = quadraticRoots(hB)[1].value, uA = x -> 0.0, uB = x -> 0.0)
    elseif hB(0) == hB(0.5) == hB(1) == 0
        return(hA = hA, hB = hB, meanA = quadraticRoots(hA)[1].value, meanB = quadraticRoots(hA)[1].value, uA = x -> 0.0, uB = x -> 0.0)
    end

    # Obtain roots of hA and hB
    rootsA = quadraticRoots(hA)
    rootsB = quadraticRoots(hB)
    T = 0.5*(hA(1)+hA(-1)-2*hA(0))
    mA = [x.multiplicity for x in rootsA]
    mB = [x.multiplicity for x in rootsB]

    # If hA and hB have distinct roots, order them so that the steady state is first
    if mA == [1,1]
        if (0 < rootsA[2].value < 1) || rootsA[1].value < 0 || rootsA[1].value > 1
            rootsA = reverse(rootsA) 
        elseif rootsA[1].value == 0 && rootsA[2].value == 1 && hA(0.5) > 0
            rootsA = reverse(rootsA)
        elseif rootsA[1].value == 1 && rootsA[2].value == 0 && hA(0.5) < 0
            rootsA = reverse(rootsA)
        end
    end
    if mB == [1,1]
        if (0 < rootsB[2].value < 1) || rootsB[1].value < 0 || rootsB[1].value > 1
            rootsB = reverse(rootsB) 
        elseif rootsB[1].value == 0 && rootsB[2].value == 1 && hB(0.5) > 0
            rootsB = reverse(rootsB)
        elseif rootsB[1].value == 1 && rootsB[2].value == 0 && hB(0.5) < 0
            rootsB = reverse(rootsB)
        end
    end

    I = sort([rootsA[1].value, rootsB[1].value])
    if isapprox(diff(I)[1], 0.0, atol=1e-10) 
        return(hA = hA, hB = hB, rootsA = rootsA, rootsB = rootsB, I = I, meanA = I[1], meanB = I[1], uA = x -> 0.0, uB = x -> 0.0)
    else
        if mA == [1]
            Aexponents = [lambda/(hA(0) - hA(1))]
        elseif mA == [1,1]
            Aexponents = [1, -1] * lambda/sqrt(quadraticDiscriminant(hA))
        elseif mA == [2]
            Aexponents = [0.0] 
        end
        if mB == [1]
            Bexponents = [lambda/(hB(0) - hB(1))]
        elseif mB == [1,1]
            Bexponents = [1, -1] * lambda/sqrt(quadraticDiscriminant(hB))
        elseif mB == [2]
            Bexponents = [0.0] 
        end
        function logQ(x)
            q = 0
            if mA == [2]
                q += -lambda/(T*(rootsA[1].value - x))
            end
            if mB == [2]
                q += -lambda/(T*(rootsB[1].value - x))
            end
            return q
        end

        # Exception for when the integral does not converge
        if mA == mB == [1,1]
            if isapprox(rootsB[2].value, rootsA[1].value, atol=1e-10)
                if -1 + Aexponents[1] + Bexponents[2] <= -1
                    return(hA = hA, hB = hB, rootsA = rootsA, rootsB = rootsB, I = I, meanA = rootsA[1].value, meanB = rootsA[1].value, uA = x -> 0.0, uB = x -> 0.0)
                end
            elseif isapprox(rootsB[1].value, rootsA[2].value, atol=1e-10)
                if -1 + Aexponents[2] + Bexponents[1] <= -1
                    return(hA = hA, hB = hB, rootsA = rootsA, rootsB = rootsB, I = I, meanA = rootsB[1].value, meanB = rootsB[1].value, uA = x -> 0.0, uB = x -> 0.0)
                end
            end
        end

        # Stationary distribution in environment A
        roots = [[r.value for r in rootsA];[r.value for r in rootsB]]
        uAexponents = [Aexponents .- mA; Bexponents]
        uAuniqueRoots = unique(roots)
        uAuniqueExponents = [sum([uAexponents[i] for i in findall(x -> isapprox(x, uAuniqueRoots[j], atol=1e-10),roots)]) for j in eachindex(uAuniqueRoots)]
        uAsingularities = findall(i -> uAuniqueRoots[i] in I && -1 < uAuniqueExponents[i] < 0, eachindex(uAuniqueRoots))   
        
        function loguA(x)
            return logQ(x) + sum([log(abs(x-uAuniqueRoots[i])) * uAuniqueExponents[i] for i in eachindex(uAuniqueRoots)])
        end
        if length(uAsingularities) == 0
            intA1 = try quadgk(x -> exp(loguA(x))*x, I[1], I[2]) catch err NaN end
            intA0 = try quadgk(x -> exp(loguA(x)), I[1], I[2]) catch err NaN end
        elseif length(uAsingularities) == 1
            function logvA(x)
                return logQ(x) + sum([log(abs(x-uAuniqueRoots[i])) * uAuniqueExponents[i] for i in eachindex(uAuniqueRoots) if i != uAsingularities[1]])
            end
            a = uAuniqueRoots[uAsingularities[1]]
            b = setdiff(I, a)[1]
            t = 1/(1+uAuniqueExponents[uAsingularities[1]])
            intA1 = try quadgk(y -> exp(logvA(a + sign(b-a)*y^t))*t*(a + sign(b-a)*y^t), 0, abs(a-b)^(1/t)) catch err NaN end
            intA0 = try quadgk(y -> exp(logvA(a + sign(b-a)*y^t))*t, 0, abs(a-b)^(1/t)) catch err NaN end
        elseif length(uAsingularities) > 1
            intA1 = NaN
            intA0 = NaN
        end

        # Stationary distribution in environment B 
        uBexponents = [Aexponents; Bexponents .- mB]
        uBuniqueRoots = unique(roots)
        uBuniqueExponents = [sum([uBexponents[i] for i in findall(x -> isapprox(x, uBuniqueRoots[j], atol=1e-10),roots)]) for j in eachindex(uBuniqueRoots)]
        uBsingularities = findall(i -> uBuniqueRoots[i] in I && -1 < uBuniqueExponents[i] < 0, eachindex(uBuniqueRoots))   

        function loguB(x)
            return logQ(x) + sum([log(abs(x-uBuniqueRoots[i])) * uBuniqueExponents[i] for i in eachindex(uBuniqueRoots)])
        end
        if length(uBsingularities) == 0
            intB1 = try quadgk(x -> exp(loguB(x))*x, I[1], I[2]) catch err NaN end
            intB0 = try quadgk(x -> exp(loguB(x)), I[1], I[2]) catch err NaN end
        elseif length(uBsingularities) == 1
            function logvB(x)
                return logQ(x) + sum([log(abs(x-uBuniqueRoots[i])) * uBuniqueExponents[i] for i in eachindex(uBuniqueRoots) if i != uBsingularities[1]])
            end
            a = uBuniqueRoots[uBsingularities[1]]
            b = setdiff(I, a)[1]
            t = 1/(1+uBuniqueExponents[uBsingularities[1]])
            intB1 = try quadgk(y -> exp(logvB(a + sign(b-a)*y^t))*t*(a + sign(b-a)*y^t), 0, abs(a-b)^(1/t)) catch err NaN end
            intB0 = try quadgk(y -> exp(logvB(a + sign(b-a)*y^t))*t, 0, abs(a-b)^(1/t)) catch err NaN end
        elseif length(uBsingularities) > 1
            intB1 = NaN
            intB0 = NaN
        end

        meanA = intA1[1]/intA0[1]
        meanB = intB1[1]/intB0[1]
        meanA = isnan(meanA) ? 0.5 : meanA
        meanB = isnan(meanB) ? 0.5 : meanB
        return(hA = hA, hB = hB, rootsA = rootsA, rootsB = rootsB, I = I, uA = x -> exp(loguA(x))/intA0[1], uB = x -> exp(loguB(x))/intB0[1], Aexponents = Aexponents, Bexponents = Bexponents, meanA = meanA, meanB = meanB)

    end
end



struct Root 
    value::Float64
    multiplicity::Int64
end

function stationaryState(mat) # stationary state of a transition matrix 
    n = nullspace(mat)
    if size(n,2) > 1
        println(mat)
        u = findfirst(i -> sameSign(n[:,i]),1:size(n,2)) # there may be multiple columns; find the first column where all entries have the same sign
        n = n[:,u]
    end
    s = n / sum(n)
    return(vec(s))
end

function quadraticDiscriminant(f) # obtain the discriminant of a quadratic
    c = f(0)
    a = 0.5*(f(1)+f(-1)-2*c)
    b = f(1)-a-c
    return b*b-4*a*c
end

function quadraticRoots(f) # obtain the roots of a quadratic, with multiplicity
    c = f(0)
    a = 0.5*(f(1)+f(-1)-2*c)
    b = f(1)-a-c
    D = b*b-4*a*c
    if isapprox(a, 0.0, atol=1e-10)
        return [Root(-c/b,1)]
    elseif isapprox(D, 0.0, atol=1e-10)
        return [Root(-b/(2*a),2)]
    else
        return [Root((-b+sqrt(D))/(2*a),1),Root((-b-sqrt(D))/(2*a),1)]
    end
end



struct Strategy3
    activityZ::Float64
    piAXY::Float64
    piAXZ::Float64
    piAYX::Float64
    piAYZ::Float64
    piAZX::Float64
    piAZY::Float64
    piBXY::Float64
    piBXZ::Float64
    piBYX::Float64
    piBYZ::Float64
    piBZX::Float64
    piBZY::Float64
    thetaXXY::Float64
    thetaXXZ::Float64
    thetaXYX::Float64
    thetaXYZ::Float64 # meet X, go from Y to Z
    thetaXZX::Float64
    thetaXZY::Float64
    thetaYXY::Float64
    thetaYXZ::Float64
    thetaYYX::Float64
    thetaYYZ::Float64
    thetaYZX::Float64
    thetaYZY::Float64
    thetaZXY::Float64
    thetaZXZ::Float64
    thetaZYX::Float64
    thetaZYZ::Float64
    thetaZZX::Float64
    thetaZZY::Float64
    rhoXY::Float64
    rhoXZ::Float64
    rhoYX::Float64
    rhoYZ::Float64
    rhoZX::Float64
    rhoZY::Float64
end

function randstruct(::Type{T}) where T
    T((rand() for FT in fieldtypes(T))...)
end

function triplets(n)
    [(x, y, n - x - y) for x in 0:n for y in 0:(n - x)]
end

storedTriplets = [triplets(n) for n in 1:100]

function tripletTransitions(n)
    t = triplets(n)
    origins = zeros(Int64, length(t),length(t))
    destinations = zeros(Int64, length(t),length(t))
    for i in 1:length(t)
        for j in 1:length(t)
            if sum((t[i] .- t[j]).^2) == 2
                origins[i,j] = findfirst(==(1), t[i] .- t[j])
                destinations[i,j] = findfirst(==(1), t[j] .- t[i])
            end
        end
    end
    return origins, destinations
end

storedTripletTransitions = [tripletTransitions(n) for n in 1:10]

function pTransition(state, originType, destinationType, env, alpha, beta, gamma, delta, s)
    p = 0.0
    states = Dict("X" => state[1], "Y" => state[2], "Z" => state[3])
    # environment response
    p += alpha * states[originType] * getproperty(s, Symbol("pi" * env * originType * destinationType))
    # interactions
    for interactionType in ["X", "Y", "Z"]
        p += beta * states[originType] * states[interactionType] * getproperty(s, Symbol("theta" * interactionType * originType * destinationType))
    end
    # random switching
    p += gamma * states[originType] * getproperty(s, Symbol("rho" * originType * destinationType))
    p += delta * states[originType]
    return p
end


function solveFiniteModel3(parameters,s)
    m = createSystemMarkovChain3(parameters["n"], parameters["alpha"], parameters["beta"], parameters["gamma"], parameters["delta"], parameters["lambda"], s)
    n = parameters["n"]
    states = triplets(n)
    k = length(states)
    activity = [Int64(round(s.activityZ)) * state[3] + state[2] for state in states]
    v = stationaryState(m)    
    uA = [(sum([activity[i] == a ? v[i] : 0 for i in 1:k])) for a in 0:n]
    uB = [(sum([activity[i] == a ? v[i+k] : 0 for i in 1:k])) for a in 0:n]
    uA = uA/sum(uA)
    uB = uB/sum(uB)
    meanA = sum(uA .* (0:n))/n
    meanB = sum(uB .* (0:n))/n 
    return(uA=uA,uB=uB,meanA=meanA,meanB=meanB)
end

function createSystemMarkovChain3(n, alpha, beta, gamma, delta, lambda, s)
    states = triplets(n)
    origins, destinations = storedTripletTransitions[n]
    k = length(states)
    m = zeros(2*k, 2*k)
    for i in 1:k
        m[i,i] += -lambda
        m[k+i, i] += lambda
        m[k+i, k+i] += -lambda
        m[i, k+i] += lambda
    end
    for i in 1:k
        for j in 1:k
            if origins[i,j] != 0
                originType = ["X", "Y", "Z"][origins[i,j]]
                destinationType = ["X", "Y", "Z"][destinations[i,j]]
                m[i,i] -= pTransition(states[i], originType, destinationType, "A", alpha, beta, gamma, delta, s)
                m[j,i] += pTransition(states[i], originType, destinationType, "A", alpha, beta, gamma, delta, s)
                m[i+k,i+k] -= pTransition(states[i], originType, destinationType, "B", alpha, beta, gamma, delta, s)
                m[j+k,i+k] += pTransition(states[i], originType, destinationType, "B", alpha, beta, gamma, delta, s)
            end
        end
    end
    return m
end

