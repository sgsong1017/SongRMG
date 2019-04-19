% Coded by Song, S. (April 2013)
% Writing an SRF source file (SRF: Standard Rupture Format by R. Graves)
%clear all

function [] = gen_srf(rup)

%scale_lat = 111.2035;  % Southern California only
%scale_lon =  92.2881;

deg2rad = pi/180.;
%%% End of inputs

disp('    ')
disp('### Step IV: Generating SRF files...')
disp('### ')

for iter = 1:rup.num
    
  outfl = [rup.name '_' num2str(iter,'%05d') '.srf'];

  %%% SRF - Data Block Input 
  for i=1:rup.nx
    for k=1:rup.nz
        
        z_km = rup.dtop + (k-0.5)*rup.dz*sin(deg2rad*rup.dip);
        
        azi = rup.stk + 90;  % fault normal direction
        len = (i - 0.5 - rup.nx/2)*rup.dx;
        
        x_km = len*sin(deg2rad*rup.stk) + ...
                   (k-0.5)*rup.dz*cos(deg2rad*rup.dip)*sin(deg2rad*azi);
        y_km = len*cos(deg2rad*rup.stk) + ...
                   (k-0.5)*rup.dz*cos(deg2rad*rup.dip)*cos(deg2rad*azi);

        rng = km2deg(sqrt(x_km^2 + y_km^2));
        azi = atan2(y_km, x_km);
        azi = (azi>=0 & azi<=pi/2).*(pi/2 - azi) + (azi > pi/2).*(2*pi+pi/2-azi) ...
                 + (azi<0).*(abs(azi) + pi/2);
        azi = rad2deg(azi);

%        dlon = x_km/scale_lon;
%        dlat = y_km/scale_lat;
        
%        subf.lon(k,i) = rup.elon + dlon;
%        subf.lat(k,i) = rup.elat + dlat;

        [subf.lat(k,i) subf.lon(k,i)] = reckon(rup.elat,rup.elon,rng,azi);

        subf.dep(k,i) = z_km;
        
        subf.stk(k,i) = rup.stk;
        subf.dip(k,i) = rup.dip;
        
        subf.area(k,i)  = (rup.dx*1e5)*(rup.dz*1e5);  % (cm * cm)
        subf.tinit(k,i) = rup.rupT.dist{iter}(k,i);
        subf.dt(k,i)    = rup.svf.dt;
        
        subf.rak(k,i)   = rup.rak;
        subf.slip1(k,i) = rup.slip.dist{iter}(k,i);
        subf.nt1(k,i)   = rup.svf.nt{iter}(k,i);
        
        % slip2 and slip3 are zeros for the moment.
        subf.slip2(k,i) = 0.;
        subf.nt2(k,i)   = 0;
        subf.slip3(k,i) = 0.;
        subf.nt3(k,i)   = 0;
        
    end
  end    
  % End of Data Block Input

  disp('   ')
  str_srf = ['=> Writing SRF file : #' num2str(iter,'%05d')];
  disp(str_srf)
  disp('   ')

  fid = fopen(outfl,'w');

  %%% writing Header Block
  fprintf(fid,'%s \n','1.0');  % Version 1.0 
  fprintf(fid,'%s \n','PLANE 1');  % Currently one segment only

  fprintf(fid,'%10.4f %10.4f %7d %7d %9.2f %9.2f \n', ...
                rup.elon,rup.elat,rup.nx,rup.nz,rup.L,rup.W);
  fprintf(fid,'%10.2f %10.2f %7.2f %7.2f %9.2f \n', ...
                rup.stk,rup.dip,rup.dtop,rup.shyp(iter),rup.dhyp(iter));
            
  %%% writing Data Block            
  fprintf(fid,'%s %7d ','POINTS',rup.nx*rup.nz);

  for k=1:rup.nz
    for i=1:rup.nx
          
      fprintf(fid,'\n %10.4f %10.4f %7.2f %10.2f %10.2f %e %7.2f %7.4f \n', ...
        subf.lon(k,i),subf.lat(k,i),subf.dep(k,i),subf.stk(k,i),subf.dip(k,i), ...
        subf.area(k,i),subf.tinit(k,i),subf.dt(k,i));
   
      fprintf(fid,'%10.2f %10.2f %6d %10.2f %6d %10.2f %6d \n', ...
        subf.rak(k,i),subf.slip1(k,i),subf.nt1(k,i), ...
        subf.slip2(k,i),subf.nt2(k,i),subf.slip3(k,i),subf.nt3(k,i));
   
      for j=1:rup.svf.nt{iter}(k,i)
        fprintf(fid,'%e',rup.svf.svf{iter}{k}{i}(j)*subf.slip1(k,i));
        fprintf(fid,'   ');
        if (mod(j,6) == 0) & (j < rup.svf.nt{iter}(k,i))
            fprintf(fid,'\n');
        end
      end  
    
%       if (mod(i,25) == 0)
%           str_i = sprintf('%10d subfaults completed',rup.nz*i);
%           disp(str_i)
%       end
    end
  end   

end
