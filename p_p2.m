%%% Coded by Song, S. (May. 2010)
%%% Analyze 2-point stats of simulated source parameters

function [rho_h rho_p] = p_p2(rup, hx, hz)

addpath ./util

dx = rup.dx; 
dz = rup.dz;

ind = {'slip', 'Vr','risT','psv'};
for iter=1:rup.num
    
    [rho_h{1}{1}{iter} rho_p{1}{1}{iter}] = co_vario2( ...
        rup.slip.dist{iter},rup.slip.dist{iter},hx,hz,[dz dx]);
    [rho_h{1}{2}{iter} rho_p{1}{2}{iter}] = co_vario2( ...
        rup.slip.dist{iter},rup.Vr.dist{iter},hx,hz,[dz dx]);
    [rho_h{1}{3}{iter} rho_p{1}{3}{iter}] = co_vario2( ...
        rup.slip.dist{iter},rup.risT.dist{iter},hx,hz,[dz dx]);
%    [rho_h{1}{4}{iter} rho_p{1}{4}{iter}] = co_vario2( ...
%       rup.slip.dist{iter},rup.psv.dist{iter},hx,hz,[dz dx]);
    
    [rho_h{2}{2}{iter} rho_p{2}{2}{iter}] = co_vario2( ...
        rup.Vr.dist{iter},rup.Vr.dist{iter},hx,hz,[dz dx]);
    [rho_h{2}{3}{iter} rho_p{2}{3}{iter}] = co_vario2( ...
        rup.Vr.dist{iter},rup.risT.dist{iter},hx,hz,[dz dx]);
%    [rho_h{2}{4}{iter} rho_p{2}{4}{iter}] = co_vario2( ...
%        rup.Vr.dist{iter},rup.psv.dist{iter},hx,hz,[dz dx]);
    
    [rho_h{3}{3}{iter} rho_p{3}{3}{iter}] = co_vario2( ...
        rup.risT.dist{iter},rup.risT.dist{iter},hx,hz,[dz dx]);
%    [rho_h{3}{4}{iter} rho_p{3}{4}{iter}] = co_vario2( ...
%        rup.risT.dist{iter},rup.psv.dist{iter},hx,hz,[dz dx]);
    
%    [rho_h{4}{4}{iter} rho_p{4}{4}{iter}] = co_vario2( ...
%        rup.psv.dist{iter},rup.psv.dist{iter},hx,hz,[dz dx]);
    
    iter
end    

for iter=1:rup.num
  for k=1:3
    for i=k:3
      rho_s{k}{i}.tmp(:,:,iter) = rho_h{k}{i}{iter};
    end
  end
end    

for k=1:3
  for i=k:3
    rho_s{k}{i}.mu = mean(rho_s{k}{i}.tmp,3);
    rho_s{k}{i}.sigma = std(rho_s{k}{i}.tmp,0,3);
  end
end  

[hxx hzz] = meshgrid(hx,hz);
for k=1:3
   for i=k:3
        
      rho_h0{k}{i} = exp(-sqrt(((hxx-rup.p2.RDx(k,i))/rup.p2.ax(k,i)).^2 ...
                             +((hzz-rup.p2.RDz(k,i))/rup.p2.az(k,i)).^2));
                         
      rho_h0{k}{i} = rup.p2.cc(k,i)*rho_h0{k}{i};
      
   end   
end        
        
%%%Plotting
cx = ceil(length(hx)/2);
cz = ceil(length(hz)/2);

figure
for iter = 1:rup.num
  for k=1:3
    for i=k:3
      subplot(4,4,(k-1)*4+i)
      plot(rho_p{k}{i}{iter}.hx,rho_h{k}{i}{iter}(cz,:),'c'), grid on, hold on
    end
  end
end  

for k=1:3
   for i=k:3
      subplot(4,4,(k-1)*4+i)
      plot(hx,rho_h0{k}{i}(cz,:),'r','linewidth',3), grid on, hold on
      
      plot(hx,rho_s{k}{i}.mu(cz,:),'b','linewidth',3)
      plot(hx,rho_s{k}{i}.mu(cz,:)-rho_s{k}{i}.sigma(cz,:),'b--','linewidth',2)
      plot(hx,rho_s{k}{i}.mu(cz,:)+rho_s{k}{i}.sigma(cz,:),'b--','linewidth',2)
      xlabel('h (km)'), ylabel('rho(h)')
      title([ind{k} ' vs. ' ind{i}])
   end
end

figure
for iter = 1:rup.num
  for k=1:3
    for i=k:3
      subplot(4,4,(k-1)*4+i)
      plot(rho_p{k}{i}{iter}.hz,rho_h{k}{i}{iter}(:,cx),'c'), grid on, hold on
    end
  end  
end
        
for k=1:3
   for i=k:3
      subplot(4,4,(k-1)*4+i)
      plot(hz,rho_h0{k}{i}(:,cx),'r','linewidth',3), grid on, hold on
      
      plot(hz,rho_s{k}{i}.mu(:,cx),'b','linewidth',3)
      plot(hz,rho_s{k}{i}.mu(:,cx)-rho_s{k}{i}.sigma(:,cx),'b--','linewidth',2)
      plot(hz,rho_s{k}{i}.mu(:,cx)+rho_s{k}{i}.sigma(:,cx),'b--','linewidth',2)
      xlabel('h (km)'), ylabel('rho(h)')
      title([ind{k} ' vs. ' ind{i}])
   end
end  



