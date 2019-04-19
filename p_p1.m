% Coded by Song, S. (April 2010)
% analyze 1-point statistics of simulated source parameters

function [slip Vr risT psv] = p_p1(rup)

slip = [];
Vr = [];
psv = [];
risT = [];

for i=1:rup.num
    slip = [slip; rup.slip.dist{i}(:)];
    Vr   = [Vr;   rup.Vr.dist{i}(:)];
    psv  = [psv;  rup.psv.dist{i}(:)];
    risT = [risT; rup.risT.dist{i}(:)];
end    

figure
subplot(221)
hist(slip,30), grid on
xlabel('slip (cm)')

subplot(222)
hist(Vr,30), grid on
xlabel('Vr (km/sec)')

subplot(223)
hist(psv,30), grid on
xlabel('Vmax (cm/sec)')

subplot(224)
hist(risT,30), grid on
xlabel('risetime (sec)')
