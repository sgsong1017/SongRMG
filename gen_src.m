% Coded by Song, S. (April. 2013)
% Basic input parameters to describe a finite source model

function [rup] = gen_src(rup)

% Basic information for an event
rup.outfl = ['rup_' rup.name '.mat'];

rup.target_Mw   = 6.6;
rup.target_Mo = fMw2MoN(rup.target_Mw);  % Nm

rup.num = 5;         % number of simulation

rup.seed.seed = 1999;
rup.seed.hypo = 1999;

rup.elon = -118.515;  % longitude, top center of the fault plane in degree
rup.elat =   34.344;  % latitude,  top center of the fault plane in degree

rup.L    = 28.20;       % fault length (along-strike) in km
rup.W    = 14.10;       % fault width (along-dip) in km

rup.stk  = 0;      % strike in deg.
rup.dip  = 90;       % dip in deg.
rup.rak  = 180;      % rake in deg.

rup.dtop = 0.;       % depth to top of fault
rand('state',rup.seed.hypo);
rup.shyp = [rand(rup.num,1) - 0.5]*rup.L*0.8;      % along strike location (from top center) of hypocenter
rup.dhyp = [rand(rup.num,1)*0.8+0.1]*rup.W;        % along dip location (from top edge) of hypocenter

rup.dx   = 0.1;      % grid size (along-strike) in km
rup.dz   = 0.1;      % grid size (along-dip) in km

rup.dx1 = 0.5;       % grid size for sub-calculation (Chol, Eig)
rup.dz1 = 0.5;       % grid size for sub-calculation

rup.svf.type = 'etinti';   % currently 'tri','rec','pliu', 'etinti' available
rup.svf.dt   = 0.02;
%%% End of inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rup.nx = ceil(rup.L/rup.dx);
rup.nz = ceil(rup.W/rup.dz);

rup.nx1 = ceil(rup.L/rup.dx1);
rup.nz1 = ceil(rup.W/rup.dz1);

for inum=1:rup.num
  rup.lx{inum} = [rup.dx/2:rup.dx:rup.dx*rup.nx] - (rup.nx*rup.dx/2 + rup.shyp(inum));
  rup.lz{inum} = [rup.dz/2:rup.dz:rup.dz*rup.nz] - rup.dhyp(inum);
  
  rup.lx1{inum} = [rup.dx1/2:rup.dx1:rup.dx1*rup.nx1] - (rup.nx1*rup.dx1/2 + rup.shyp(inum));
  rup.lz1{inum} = [rup.dz1/2:rup.dz1:rup.dz1*rup.nz1] - rup.dhyp(inum);
  
  rup.lx1{inum} = [(rup.lx1{inum}(1)-rup.dx1) rup.lx1{inum} (rup.lx1{inum}(end)+rup.dx1)];
  rup.lz1{inum} = [(rup.lz1{inum}(1)-rup.dz1) rup.lz1{inum} (rup.lz1{inum}(end)+rup.dz1)];
  
  [XX ZZ] = meshgrid(rup.lx{inum},rup.lz{inum});
  rup.dis{inum} = sqrt(XX.^2 + ZZ.^2);  % distance from the nucleation point on the fault plane

end

  rup.nx1 = rup.nx1 + 2;
  rup.nz1 = rup.nz1 + 2;
  


