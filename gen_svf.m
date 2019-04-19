% Coded by Song, S. (December 2009)
% Generating triangular or rectangular shape Slip Velocity Function
% Input : slip duration(tw), dt, nt, svf type = ['tri' or 'rec']
% Output: svf, time

function [svf, time] = gen_svf(tw,dt,nt,type,par)

switch type

    case 'tri'
    
        h = 2/tw;
        m = h/(tw/2);
      
        time = [0:dt:dt*(nt-1)];
        
        svf1 = m*time(time <= tw/2); 
        svf2 = m*tw/2 - m*(time(time > tw/2 & time <= tw) - tw/2);
      
        svf = zeros(1,nt);
        svf3 = [svf1 svf2];
        
        if (time(end) < tw)
            
            disp('not enough time steps')
            svf = nan;
            
        else
            
            svf(1:length(svf3)) = svf3;
            
        end    

    case 'rec'
      
        h = 1/tw;
        time = [0:dt:dt*(nt-1)];
     
        svf  = zeros(1,nt);
        
        if (time(end) < tw)
            
            disp('not enough time steps');
            svf = nan;
         
        else
            
            svf(time < tw) = h;
            
        end    
    
    case 'pliu'   % based on svf by Liu et al. (BSSA, 2006)
        
        [svf, time] = svf_pliu(tw,dt,nt);
        
    case 'etinti' % based on svf by Tinti et al. (BSSA, 2005)

        T_acc = par.T_acc;  % fixed T_acc
%        T_acc = tw*0.1;      % 10% of total risetime

        tau_s = T_acc/1.3;
        tau_r = tw - 2*tau_s;

        if (tau_r <= tau_s)

           tau_r = tau_s*1.01;

        end
 
        [svf, time] = svf_etinti(tau_s,tau_r,dt,nt);

    otherwise
    
        disp('Wrong SVF type');
        
end        
      
