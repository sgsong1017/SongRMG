%%% Coded by Song, S. (Feb. 2010)
%%% Generate spatial coherence structure...
%%% Currently exponential function only

function [rho] = gen_rho(rup,r)

N = length(rup.p1.mu);

for k=1:N
    for i=k:N
        
        rho{k}{i} = exp(-sqrt(((r.x-rup.p2.RDx(k,i))/rup.p2.ax(k,i)).^2 ...
                             +((r.z-rup.p2.RDz(k,i))/rup.p2.az(k,i)).^2));
                         
        rho{k}{i} = rup.p2.cc(k,i)*rho{k}{i};
        
    end
end    

  
