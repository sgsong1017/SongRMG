%%% Coded by Song, S. (April 2013)
%%% Generate 2D distributions of source parameters
%%% using the Cholesky factorization
%%% based on 1-point and 2-point statistics

function [rup] = gen_dist(rup)

tic

N  = rup.nx*rup.nz;
N1 = (rup.nx1)*(rup.nz1);

str_N = sprintf('=> Number of subfaults (fine grid): %d',N);
disp(str_N), disp('  ')

str_N = sprintf('=> Number of subfaults (coarse grid): %d',N1);
disp(str_N), disp('  ')

[XX ZZ] = meshgrid(rup.lx1{1},rup.lz1{1});

X = single(XX(:));
Z = single(ZZ(:));

Z1 = repmat(Z ,1,N1); clear Z
r.z = (Z1' - Z1);    clear Z1  
        
X1 = repmat(X ,1,N1); clear X
r.x = (X1' - X1);    clear X1  

rho = gen_rho(rup,r); %clear r

Cm = [rho{1}{1}  rho{1}{2}  rho{1}{3}; ...
      rho{1}{2}' rho{2}{2}  rho{2}{3}; ...
      rho{1}{3}' rho{2}{3}' rho{3}{3}];  %clear rho

t1 = toc; 
str_t1 = sprintf('=> Elapsed time for constructing Cm: %10.2f',t1);
disp(str_t1)

%removing negative eigen values for positivity constraints
[eig_v eig_d] = eig(Cm);
rup.eigen = diag(eig_d);

eigen = rup.eigen;
eigen(eigen<0.01) = 0.01;

Cm = eig_v*diag(eigen)*(eig_v)';

t2 = toc; 
str_t2 = sprintf('=> Elapsed time for Eigen decomposition: %10.2f',t2-t1);
disp(str_t2)


L = chol(Cm,'lower');   clear Cm

t3 = toc; 
str_t3 = sprintf('=> Elapsed time for Cholesky Factorization: %10.2f',t3-t2);
disp(str_t3)

randn('state',rup.seed.seed);
%randn('state',sum(100*clock));
for iter=1:rup.num
  s0 = randn(3*N1,1);
  s = L*s0; %clear L

  slip1 = s(1:N1);
  Vr1   = s(N1+1:2*N1);
  psv1 = s(2*N1+1:end);

  slip1 = reshape(slip1,rup.nz1,rup.nx1);
  Vr1   = reshape(Vr1  ,rup.nz1,rup.nx1);
  psv1  = reshape(psv1,rup.nz1,rup.nx1);

  rup.slip1.dist{iter} = slip1;  
  rup.Vr1.dist{iter}   = Vr1;    
  rup.psv1.dist{iter} = psv1;
end  

t4 = toc; 
str_t4 = sprintf('=> Elapsed time for random sampling: %10.2f',t4-t3);
disp(str_t4)


