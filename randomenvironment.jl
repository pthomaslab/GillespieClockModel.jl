### Simulation of cyanobacterial circadian clock dynamics in time-varying enviroments
## mass-action version adapted from Chew et al. Nature Communications 9: 3004 (2018)
## here it the phosphorylation rate is assumed to be modulated by light levels 

using DataFrames, CSV, Plots

# load the simulator
include("sim.jl")

# some parameters
NRepeats = 20

# computes phosphorylation state
function getPhos(df; offset=2)
    pmax = 5
    pd = reduce(hcat,[sum(eachcol(df[:, offset + 4*p:offset + 4*p+3])) for p in 0:pmax])
    pp = reduce(hcat,[(p/pmax)*sum(eachcol(df[:, offset + 4*p:offset + 4*p+3])) for p in 0:pmax])
    gg = vec(sum(pp ./ sum(pd,dims=2),dims=2))
    time = df[:,1]

    time, gg
end

# maps instanteous light levels to phosphorylation rate via Michaelis-Menten type relation (fitted from medium and low LD data)
function MM(x;vm=.6,K=0.8,c=0.05)
    ldVal=vm # set level 
    v=-(((c - ldVal)*(K + ldVal))/ldVal)
    return (v*x/(x+K)+c)
end

function makePlot(fnameTime,fnameLight,ptitle)
    # load the model
    jprob = createJumpProblem()

    # import environment
    tt=Matrix(CSV.read(fnameTime,DataFrame,header=false))[1,:]
    ld=Matrix(CSV.read(fnameLight,DataFrame,header=false))[1,:]
    env0=reduce(hcat,[tt,ld])
    light0=reduce(hcat,[tt,ld])

    # apply mapping to phosphorylation rate
    env0[:,2] .= MM.(env0[:,2])

    # simulate and plot
    scale = 0.3 /maximum(light0[:,2]) # scaling factor for visualisation
    p1=plot(light0[:,1],scale  * light0[:,2], color=:grey, fill=(0, 0.9, :gold),xlabel="time [h]",ylabel="phosphorylation level", title=ptitle)
    for n in 1:NRepeats
        df =simulateLineage(env0,jprob)
        time, phos = getPhos(df)
        plot!(p1,time,phos; color=:deepskyblue, legend=false, xlimits=(0,maximum(light0[:,1])))
    end
    return p1
end

### Set working directory

cd("/home/pthomas/Dropbox/clockModel/test_LLfit/final/")

### 1) simulate the noisy day start and end environments

p1 = makePlot("lightinputs/noisydaystart_time.csv","lightinputs/noisydaystart_light.csv","noisy day start")
p2 = makePlot("lightinputs/noisydayend_time.csv","lightinputs/noisydayend_light.csv","noisy day end")

### 2) simulate the natural environment

p3 = makePlot("lightinputs/c1_time.csv","lightinputs/c1_light0.csv","Caribbean 1")
p4 = makePlot("lightinputs/c2_time.csv","lightinputs/c2_light.csv","Caribbean 2")

### make combined plot
plot(p1,p2,p3,p4,layout=[[2],[2]])
savefig("out/output.pdf")