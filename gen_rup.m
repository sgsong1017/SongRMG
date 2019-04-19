% Coded by Song, S. (April 2013)
% Rupture model generator based on 1-point and 2-point statistics 
% of key kinematic source parameters
%clear all

function [rup] = gen_rup(rup)

deg2rad = pi/180.;

par.tag = 'Additional Parameters';

addpath ./FastMarching/code
addpath ./FastMarching/code/functions

% Step I: 2-point statistics
% generate 2D distributions of source parameters 
% assuming multi-variate Gaussian with zero mean and unit std
disp('   ')
disp('### Step I: 2-point statistics   ')
disp('### generating 2D distributions of source parameters')
disp('### based on covariance matrix constructed from     ')
disp('### auto- and cross-correlation ')
disp('   ')

rup = gen_dist(rup);

%%%Interpolation
for iter=1:rup.num
   
    [X Z] = meshgrid(rup.lx{iter},rup.lz{iter});
    [X1 Z1] = meshgrid(rup.lx1{iter},rup.lz1{iter});
    
    rup.slip.dist{iter} = interp2(X1,Z1,rup.slip1.dist{iter},X,Z);
    rup.Vr.dist{iter}   = interp2(X1,Z1,rup.Vr1.dist{iter},X,Z);
    rup.psv.dist{iter}  = interp2(X1,Z1,rup.psv1.dist{iter},X,Z);
    
end

% Step II: 1-point statistics
% currently we assume the Gaussian distribution 
disp('   ')
disp('### Step II: 1-point statistics   ')
disp('### Adjust 1-point statistics, e.g., mean and sigma')
disp('### Currently assume Gaussian distribution         ')
disp('### But it can be transformed into non Gaussian distribution')
disp('### if necessary   ')
disp('   ')
tic

  
for iter=1:rup.num
  
  %%% Slip  
  %fixing mean slip and Mo    
  rup.slip.dist{iter} = (rup.slip.dist{iter}-mean2(rup.slip.dist{iter}))/std2(rup.slip.dist{iter});  
  rup.slip.dist{iter} = rup.p1.sig(1)*rup.slip.dist{iter} + rup.p1.mu(1); 
  
  %boundary tapering for slip
  sfac = b_taper(rup.nz,rup.nx,rup.dtop);
  rup.slip.dist{iter} = rup.slip.dist{iter}.*sfac;
  
  rup.slip.dist{iter} = (rup.slip.dist{iter}-mean2(rup.slip.dist{iter}))/std2(rup.slip.dist{iter});  
  rup.slip.dist{iter} = rup.p1.sig(1)*rup.slip.dist{iter} + rup.p1.mu(1); 
  rup.slip.dist{iter}(rup.slip.dist{iter}<rup.p1.min(1)) = rup.p1.min(1);
  rup.slip.dist{iter}(rup.slip.dist{iter}>rup.p1.max(1)) = rup.p1.max(1);
  rup.p1.mu1(iter,1) = mean2(rup.slip.dist{iter});
  rup.p1.sig1(iter,1) = std2(rup.slip.dist{iter});

  %%% Rupture Velocity
  rup.Vr.dist{iter} = rup.p1.sig(2)*rup.Vr.dist{iter} + rup.p1.mu(2);
  rup.Vr.dist{iter}(rup.Vr.dist{iter}<rup.p1.min(2)) = rup.p1.min(2);
  rup.Vr.dist{iter}(rup.Vr.dist{iter}>rup.p1.max(2)) = rup.p1.max(2);
  rup.p1.mu1(iter,2) = mean2(rup.Vr.dist{iter});
  rup.p1.sig1(iter,2) = std2(rup.Vr.dist{iter});

  %%% Peak Slip Velocity
  rup.psv.dist{iter} = rup.p1.sig(3)*rup.psv.dist{iter} + rup.p1.mu(3);
  rup.psv.dist{iter}(rup.psv.dist{iter}<rup.p1.min(3)) = rup.p1.min(3);
  rup.psv.dist{iter}(rup.psv.dist{iter}>rup.p1.max(3)) = rup.p1.max(3);
  rup.p1.mu1(iter,3) = mean2(rup.psv.dist{iter});
  rup.p1.sig1(iter,3) = std2(rup.psv.dist{iter});
 
  switch rup.svf.type

  case 'rec'
  
    rup.risT.dist{iter} = rup.slip.dist{iter}./rup.psv.dist{iter};
    rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
    rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
    rup.p1.mu1(iter,4) = mean2(rup.risT.dist{iter});
    rup.p1.sig1(iter,4) = std2(rup.risT.dist{iter});
    
  case 'tri'
    
    rup.risT.dist{iter} = rup.slip.dist{iter}./rup.psv.dist{iter}*2;
    rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
    rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
    rup.p1.mu1(iter,4) = mean2(rup.risT.dist{iter});
    rup.p1.sig1(iter,4) = std2(rup.risT.dist{iter});        
    
  case 'pliu'

    C_pliu = 1.4*pi*0.13 + 1.2*0.13 + 0.3*pi*0.87;
    rup.risT.dist{iter} = 2*rup.slip.dist{iter}*pi./rup.psv.dist{iter}/C_pliu;
    rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
    rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
    rup.p1.mu1(iter,4) = mean2(rup.risT.dist{iter});
    rup.p1.sig1(iter,4) = std2(rup.risT.dist{iter});

  case 'etinti'
    
    par.T_acc = 0.2;
    rup.risT.dist{iter} = (1.04*rup.slip.dist{iter}./(par.T_acc^0.54*rup.psv.dist{iter})).^(1/0.47);  
    rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
    rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
    rup.p1.mu1(iter,4) = mean2(rup.risT.dist{iter});
    rup.p1.sig1(iter,4) = std2(rup.risT.dist{iter});

    rup.svf.T_acc = par.T_acc;  %check T_acc inconsistency (09/2015)

  otherwise

    error('The type of svf is not supported yet')

  end

  [rup.Mo(iter) rup.Mw(iter)] = fmomentN(rup.slip.dist{iter},[rup.nz*rup.dz rup.nx*rup.dx]);
  
end
t1 = toc; 
str_t1 = sprintf('=> Elapsed time for 1-point adjustment: %10.2f',t1);
disp(str_t1);


% Step III: post-processing of generated source models
% a)rupture time calculation from simulated local rupture velocity
disp('   ')
disp('### Step III: Post-processing   ')
disp('###   ')
for iter=1:rup.num
%    rup.rupT.dist{iter} = gen_rupT(rup.Vr.dist{iter}, ...
%                                rup.lx{iter},rup.lz{iter});
    shypo = rup.nx*rup.dx/2 + rup.shyp(iter);
    dhypo = rup.dhyp(iter);
    rup.rupT.dist{iter} = fm(rup.Vr.dist{iter},[shypo dhypo]',[rup.dx rup.dz], ...
        struct('implementation','MATLAB','useMC',0,'order',2));
end
t2 = toc; 
str_t2 = sprintf('=> Elapsed time for computing rupture time distribution: %10.2f',t2-t1);
disp(str_t2);


% b) generating slip velocity functions (SVF)
for iter=1:rup.num
  for k=1:rup.nz
    for i=1:rup.nx
      rup.svf.nt{iter}(k,i) = ceil((rup.risT.dist{iter}(k,i)+0.5)/rup.svf.dt); 
      [rup.svf.svf{iter}{k}{i} rup.svf.time{iter}{k}{i}] = ...
        gen_svf(rup.risT.dist{iter}(k,i),rup.svf.dt,rup.svf.nt{iter}(k,i),rup.svf.type,par);
    
      rup.psv1.dist{iter}(k,i) = max(rup.svf.svf{iter}{k}{i})*rup.slip.dist{iter}(k,i);
    end
  end
end    


t3 = toc; 
str_t3 = sprintf('=> Elapsed time for generating slip velocity functions: %10.2f',t3-t2);
disp(str_t3);


% c) generating moment rate function
Area = rup.dx*rup.dz*1e6;
for iter=1:rup.num
  tlen = ceil(max(rup.rupT.dist{iter}(:))+max(rup.risT.dist{iter}(:)))+1;
  nt = ceil(tlen/rup.svf.dt);
  mrf_time = [0:rup.svf.dt:rup.svf.dt*(nt-1)];
  mrf = zeros(rup.nx*rup.nz,nt);

  ind = 0;
  for k=1:rup.nz
    for i=1:rup.nx
      ind = ind + 1;
      i1 = round(rup.rupT.dist{iter}(k,i)/rup.svf.dt)+1;
      i2 = i1 + length(rup.svf.svf{iter}{k}{i}) - 1;
      mrf(ind,i1:i2) = rup.mu*Area*rup.svf.svf{iter}{k}{i}*rup.slip.dist{iter}(k,i)/100.;
    end
  end
  rup.mrf.mrf{iter}  = sum(mrf);
  rup.mrf.time{iter} = mrf_time;
end

t4 = toc;
str_t4 = sprintf('=> Elapsed time for generating moment rate function: %10.2f',t4-t3);
disp(str_t4);

%%% Save generated rupture models
save(rup.outfl,'rup','-v7.3');





