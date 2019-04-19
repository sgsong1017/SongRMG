% Coded by Song, S. (April 2010)
% Generating the shape of slip velocity function (SVF)
% based on Liu et al. (BSSA, 2006)

function [svf, t] = svf_pliu(tau,dt,nt)

% tau = 10;  % slip duration or rise time

t = [0:dt:dt*(nt-1)];

tau1 = 0.13*tau;
tau2 = tau - tau1;
%%% End of inputs


Cn = pi/(1.4*pi*tau1 + 1.2*tau1 + 0.3*pi*tau2);

svf = Cn*(0.7 - 0.7*cos(pi*t/tau1) + 0.6*sin(0.5*pi*t/tau1)).*(t >= 0 & t < tau1) ...
    + Cn*(1.0 - 0.7*cos(pi*t/tau1) + 0.3*cos(pi*(t - tau1)/tau2)).*(t >= tau1 & t < 2*tau1) ...
    + Cn*(0.3 + 0.3*cos(pi*(t-tau1)/tau2)).*(t >= 2*tau1 & t < tau);
  
 
%plot(t,svf,'o-'), grid on  
