% Coded by Song, S. (August 2012)
% Generating the shape of slip velocity function (SVF)
% based on Tinti et al. (BSSA, 2005)

function [svf, t] = svf_etinti(tau_s,tau_r,dt,nt)


t = [0:dt:dt*(nt-1)];

if (tau_r <= tau_s)
    
    disp('Tau_r should be larger than Tau_s');
    return;
    
end    

K = 2/(pi*tau_r*tau_s^2);  % U_tot (final slip) is one here

C1 = (1/2*t + 1/4*tau_r).*sqrt(t.*(tau_r-t)) + (t*tau_r - tau_r^2) ...
    .*asin(sqrt(t/tau_r)) - 3/4*tau_r^2*atan(sqrt((tau_r - t)./t));

C2 = 3/8*pi*tau_r^2;

C3 = (tau_s - t - 1/2*tau_r).*sqrt((t - tau_s).*(tau_r - t + tau_s)) ...
    + tau_r*(2*tau_r - 2*t + 2*tau_s).*asin(sqrt((t-tau_s)/tau_r)) ...
    + 3/2*tau_r^2*atan(sqrt((tau_r - t + tau_s)./(t - tau_s)));

C4 = (-tau_s + 1/2*t + 1/4*tau_r).*sqrt((t - 2*tau_s).*(tau_r - t + 2*tau_s)) ...
    + tau_r*(-tau_r + t - 2*tau_s).*asin(sqrt((t - 2*tau_s)/tau_r)) ...
    - 3/4*tau_r^2*atan(sqrt((tau_r - t + 2*tau_s)./(t - 2*tau_s)));

C5 = pi/2*tau_r*(t-tau_r);

C6 = pi/2*tau_r*(2*tau_s - t + tau_r);


if (tau_r > 2*tau_s )
    
    svf = (C1 + C2).*(t>=0 & t<tau_s) ...
        + (C1 - C2 + C3).*(t>=tau_s & t<2*tau_s) ...
        + (C1 + C3 + C4).*(t>=2*tau_s & t<tau_r) ...
        + (C5 + C3 + C4).*(t>=tau_r & t<tau_r+tau_s) ...
        + (C4 + C6).*(t>=(tau_r+tau_s) & t<(tau_r+2*tau_s)) ...
        + 0.*(t>=(tau_r+2*tau_s));
    
    svf = K*svf;

else
    
   svf =  (C1 + C2).*(t>=0 & t<tau_s) ...
        + (C1 - C2 + C3).*(t>=tau_s & t<tau_r) ...
        + (C5 + C3 - C2).*(t>=tau_r & t<2*tau_s) ...
        + (C5 + C3 + C4).*(t>=2*tau_s & t<tau_r+tau_s) ...
        + (C4 + C6).*(t>=(tau_r+tau_s) & t<(tau_r+2*tau_s)) ...
        + 0.*(t>=(tau_r+2*tau_s)); 
    
   svf = K*svf;
   
end    


