% Coded by Song, S.G. (August 2009)
% Computing 2D correlogram, rho(h), from 2D distributions
% Treat mean and variance of head and tail points separately
% based on P. Goovaerts (1997, p28).
%
%
% Input  : 
%           V1 : Tail distribution
%           V2 : Head distribution
%           hx : horizontial estimation range 
%           hz : vertical estimation range
%           samp = [dz dx] : grid spacing
% Output :
%           rho_h : correlogram estimates
%           p     : some useful output information
%             p.hx    : horizontal estimation range
%             p.hz    : vertical estimation range
%             p.ns(k,i) : number of samples used for each estimation
%                             rho_h(k,i)
%             p.m1 : mean of V1
%             p.m2 : mean of V2
%             p.sigma1 : std of V1
%             p.sigma2 : std of V2
%             p.ms1, p.ms2, p.sigs1, p.sigs2, separation vector h specific

function [rho_h, p] = co_vario2(V1,V2,hx,hz,samp)

[nz nx] = size(V1);

dz = samp(1); 
dx = samp(2); 

% dhz = hz(2) - hz(1);
% dhx = hx(2) - hx(1);
% 
% if (dhz < dz)
%     hz = [hz(1):dz:hz(end)];
% end
% 
% if (dhx < dx)
%     hx = [hx(1):dx:hx(end)];
% end    
    
Nhz = length(hz);
Nhx = length(hx); 

tmp1 = V1; 
tmp2 = V2;

tmp1(isnan(V1) | isnan(V2)) = nan;
tmp2(isnan(V1) | isnan(V2)) = nan;

V1 = tmp1; V2 = tmp2;

sigma1 = nanstd(V1(:));
sigma2 = nanstd(V2(:));

C_0 = sigma1*sigma2;

m1 = nanmean(V1(:)); 
m2 = nanmean(V2(:));

for jx=1:Nhx
  for jz = 1:Nhz
      
      for i=1:nx
        for k=1:nz
            
           ihx = i+round(hx(jx)/dx);
           khz = k+round(hz(jz)/dz);
        
           if (ihx >= 1 & ihx <= nx) & (khz >= 1 & khz <= nz)
               
              tmp(k,i) = V1(k,i)*V2(khz,ihx);
              t_V1(k,i) = V1(k,i);
              t_V2(k,i) = V2(khz,ihx);
              
           else
               
              tmp(k,i) = nan;
              t_V1(k,i) = nan;
              t_V2(k,i) = nan;
              
           end    
           
        end
      end
      
      ind = ~isnan(tmp(:));
      
      p.ns(jz,jx) = sum(ind);
      
      p.sample1{jz}{jx} = t_V1(ind);
      p.sample2{jz}{jx} = t_V2(ind);
      
      p.ms1(jz,jx) = mean(t_V1(ind));
      p.ms2(jz,jx) = mean(t_V2(ind));

      p.sigs1(jz,jx) = std(t_V1(ind));
      p.sigs2(jz,jx) = std(t_V2(ind));

      C_h(jz,jx) = mean(tmp(ind)); 
      
  end    
end    

C_h = C_h - p.ms1.*p.ms2;
rho_h = C_h./p.sigs1./p.sigs2;    % Covariance function to Correlogram

p.hx = hx;
p.hz = hz;

p.m1 = m1;
p.m2 = m2;

p.sigma1 = sigma1;
p.sigma2 = sigma2;

