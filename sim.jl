### Simulation of circadian clock model in time-dependent enviroment

using JumpProcesses

# import the model definition
include("clockmodel.jl")

# preset light level for entrainment
ldVal = 1.2

# initialise Gillespie simulation
function createJumpProblem(;m=5)

    rs = produceMAModel()

    Scaling=1200. # parameter setting the KaiC copy number
    sc = 1. # light scaling
    N=m+1 # number of phosphorylation states
    p = [0.429,0.49,3.57E-4*sc,2.14E-2,1.62E-4] # rate constants

    # time interval
    tspan = (0.,1000*24.)

    # initial conditions
    u0 = [0. for i in 1:4*N+1]
    u0[1]=2.0*Scaling
    u0[end]=2.0*Scaling

    println(p)

    # stochastic simulation
    dprob = DiscreteProblem(rs, u0, tspan, p)
    return JumpProblem(rs, dprob, Direct(); save_positions = (false, false))

end

# run simulator for constant light level (sc)
function runLD(solIn,jprob,tspan; sc=1.)

    p = [0.429,0.49*sc,3.57E-4,2.14E-2,1.62E-4] # rate constants

    # initial condition
    u0 = solIn[end]

    jprob2 = remake(jprob, u0 = u0, p = p, tspan = tspan)
    sol = solve(jprob2, SSAStepper())
    
    append!(solIn,[sol[end]])
end


# entrain simulator
function entrain(jprob;m=5,Scaling=1200.,ldVal=ldVal)
    N=m+1
    u0 = [0. for i in 1:4*N+1]
    u0[1]=2.0*Scaling
    u0[end]=2.0*Scaling
    sol = [u0]
    # entrain using LD square waves
    for day in 0:10
        for dd in [1,0] #day/night
            ld = (dd>0. ? ldVal : 0.0)
            tspan2 = ((day+dd)*24,(day+dd+0.5)*24)
            runLD(sol,jprob,tspan2; sc=MM.(ld))
        end 
    end
    return sol
end

# entrain simulator and return entrained initial condition
function makeICs(jprob;m=5,Scaling=1200.,ldVal=ldVal)
    sol=entrain(jprob;m=m,Scaling=Scaling,ldVal=ldVal)
    return [sol[end]]
end

# simulate lineage in a given environment (env)
function simulateLineage(env,jprob;ldVal=ldVal)

    # entrain first
    sol = makeICs(jprob,ldVal=ldVal)
    # now run the enviroment
    time = [] 
    append!(time,0.)
    step=1
    for idx in 1:step:(size(env,1)-step)
        if env[idx,1] <= 0
            continue
        end
        tspan2=(env[idx,1],env[idx+step,1])
        ldval = env[idx+step,2] 
        append!(time,env[idx,1])
        runLD(sol,jprob,tspan2; sc=ldval)
    end

    for idx in 1:length(sol)
        prepend!(sol[idx],time[idx])
    end

    df = DataFrame(transpose(reduce(hcat,sol)),:auto)

    return df
end
