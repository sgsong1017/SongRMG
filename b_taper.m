%%% Coded by Songv, S. (May 2014)
% Tapering slip values near the boundary

function [sfac] = b_taper(Nz,Nx,dtop)

%Nz = 199;
%Nx = 300;
%dtop = 0;

nz_bt = round(Nz/5);
nx_bt = round(Nx/5);

sfac = ones(Nz,Nx);
  
x = linspace(0,1,nx_bt);
z = linspace(0,1,nz_bt)';

sfac_x = sqrt(1-x.^2);
sfac_z = sqrt(1-z.^2);


if dtop < 0.1
    
  sfac_x = [repmat(fliplr(sfac_x),Nz,1) ones(Nz,Nx-2*nx_bt) repmat(sfac_x,Nz,1)];
  sfac_z = [ones(Nz-nz_bt,Nx); repmat(sfac_z,1,Nx)]; 
    
else
   
  sfac_x = [repmat(fliplr(sfac_x),Nz,1) ones(Nz,Nx-2*nx_bt) repmat(sfac_x,Nz,1)];
  sfac_z = [repmat(flipud(sfac_z),1,Nx); ones(Nz-2*nz_bt,Nx); repmat(sfac_z,1,Nx)];
  
end

sfac = sfac_x.*sfac_z;

%imagesc(sfac), axis equal tight, colorbar
    
    
    
    
    
