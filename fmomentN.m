function [Mo,Mw] = fmomentN(in,dim)
%
% function [Mo,Mw] = fmomentN(in,dim) calculates the 
% moment of an event for a given slip distribution (in cm)
% if only IN is given, DIM = size(IN)
%
% INPUT:  in 	- 2D-array that constains slip values
%	      dim	- source dimensions [W L]
%
% OUTPUT: Mo	- seismic moment, in Nm
%	      Mw	- moment magnitude
%
% mmai, 02/08/98
% --------------

if nargin == 1
  dim = size(in); 
end
W = dim(1); L = dim(2);

mu = 3.3*1e10;			% shear modulus = rigidity [N/m^2]
s = mean(in(:))/100;		% average slip over fault surface [m]
A = L*W*1e6;			% fault surface [m]

Mo = mu*s*A;			% moment [Nm]
Mw = fMo2MwN(Mo);		% moment magnitude
