function [Mw] = fMo2MwN(Mo)
%
%	function [Mw] = fMo2MwN(Mo) calculates the
%	magnitude Mw for a given moment in Nm
%	
%	for reference, see Thorne & Lay, page 384
%
%	mmai, 02/08/98
%	--------------

Mw = (2/3)*(log10(Mo) - 9.05);
