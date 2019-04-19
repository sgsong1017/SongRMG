% Coded by Song, S. (April 2010)
% plotting simulated rupture models 

function [] = p_rup(rup,iter)

figure
subplot(511)
imagesc(rup.lx{iter},rup.lz{iter},rup.slip.dist{iter}), axis equal tight, colorbar
title('slip (cm)')

subplot(512)
T_inc = 1;
[c, h] = contour(rup.lx{iter},rup.lz{iter},rup.rupT.dist{iter},[T_inc:T_inc:max(rup.rupT.dist{iter}(:))]); 
axis equal tight ij, colorbar, clabel(c,h);
title('rupture time (sec)')

subplot(513)
imagesc(rup.lx{iter},rup.lz{iter},rup.Vr.dist{iter}), axis equal tight, colorbar
title('rupture velocity (km/sec)')

subplot(514)
imagesc(rup.lx{iter},rup.lz{iter},rup.psv.dist{iter}), axis equal tight, colorbar
title('Vmax (cm/sec)')

subplot(515)
imagesc(rup.lx{iter},rup.lz{iter},rup.risT.dist{iter}), axis equal tight, colorbar
title('slip duration (sec)')
xlabel('along-strike (km)'), ylabel('down-dip (km)')
