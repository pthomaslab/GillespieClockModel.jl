using ModelingToolkit, Catalyst

### stochastic model of the core phosphorylation loop of the cyanobacterial circadian clock 
### mass-action version adapted from Chew et al. Nature Communications 9: 3004 (2018)
### here it the phosphorylation rate is assumed to be modulated by light levels 

  function produceMAModel(;m=5)
  
  # number of phosphorylation states
  N=m+1
  
  # parameters, species
  @parameters kdephos kphos kAon kAoff kABC
  @variables t
  @species (C(t))[1:N] (AC(t))[1:N] (BC(t))[1:N] (ABC(t))[1:N] A(t)
  
  rs = make_empty_network()
  
  addparam!(rs,kdephos)
  addparam!(rs,kphos)
  addparam!(rs,kAon)
  addparam!(rs,kAoff)
  addparam!(rs,kABC)
  
  for i in 1:N
      addspecies!(rs,C[i])
      addspecies!(rs,AC[i])
      addspecies!(rs,BC[i])
      addspecies!(rs,ABC[i])
  end
  addspecies!(rs,A)
  
  ## define reactions

  # spontaneous dephosphorylation
  for i in 2:N
      addreaction!(rs,Reaction(kdephos, [C[i]], [C[i-1]],[1],[1]))
  end
  # KaiA-dependent phosphorylation
  for i in 1:N-2
      addreaction!(rs,Reaction(kphos, [AC[i]], [AC[i+1]],[1],[1]))
  end
  addreaction!(rs,Reaction(kphos, [AC[N-1]], [BC[N],A],[1],[1,1]))
  # KaiA binding/unbinding
  for i in 1:N
      addreaction!(rs,Reaction(kAon, [C[i],A], [AC[i]],[1,1],[1]))
      addreaction!(rs,Reaction(kAoff, [AC[i]], [C[i],A],[1],[1,1]))
      addreaction!(rs,Reaction(kABC, [BC[i],A], [ABC[i]],[1,6],[1]))
  end
  # KaiB-dependent dephosphorylation
  for i in 3:N
    addreaction!(rs,Reaction(kdephos, [BC[i]], [BC[i-1]],[1],[1]))
    addreaction!(rs,Reaction(kdephos, [ABC[i]], [ABC[i-1]],[1],[1]))
  end
  begin
      addreaction!(rs,Reaction(kdephos, [BC[N]], [C[N-1]],[1],[1]))
      addreaction!(rs,Reaction(kdephos, [ABC[2]], [C[1],A],[1],[1,6]))
      addreaction!(rs,Reaction(kdephos, [BC[2]], [C[1]],[1],[1]))
  end
  
  return rs
  
  end