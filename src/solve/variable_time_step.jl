function variable_time_step(periods::Vector,Nsteps::Vector)
    np = size(periods,1)
    ns = size(Nsteps,1)
    ns != np &&  (error(" periods and Nsteps must have equal size"))
    Npas = sum(Nsteps) 
    j = 1
    Times = zeros(Npas+1)  
    h = zeros(Npas)
    for i = 1:ns 
        hi = periods[i]/Nsteps[i] 
        for k = 1:Nsteps[i]
            h[j] = hi
            j += 1
            Times[j] = Times[j-1] + hi
            
        end
    end
    Times[Npas+1] = Times[Npas] + periods[ns]/Nsteps[ns] 
    
    return Times, h
end 