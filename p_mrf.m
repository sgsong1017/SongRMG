%%% Coded by Song, S.G. (March 2016)
%%% Plotting Moment rate function and its spectrum...


function [nt, dt] = p_mrf(rup, irup)

dt = rup.svf.dt;
nt = length(rup.mrf.time{irup});

amp_mrf = abs(fft(rup.mrf.mrf{irup}));

amp_time = [0:1/(nt*dt):1/(nt*dt)*(nt-1)];


%%% Plotting....
figure
subplot(211)
plot(rup.mrf.time{irup},rup.mrf.mrf{irup},'b','linewidth',2), grid on
title('Moment rate function')
xlabel('time (sec)'), ylabel('M(t), Nm/sec')
set(gca,'fontsize',14);

subplot(212)
loglog(amp_time,amp_mrf,'b','linewidth',2), grid on, hold on
xlim([0 1/2/dt])
title('Amplitude spectrum')
xlabel('Hz'), ylabel('Amplitude')

loglog([0.1 10],[1e20 1e16],'r','linewidth',3)

set(gca,'fontsize',14);

