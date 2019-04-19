function [Mo] = fMw2MoN(Mw)
%
%	function [Mo] = fMw2MoN(Mw) calculates the
%	moment for a given magnitude [in Nm]
%
%	for reference see Thorne & Lay, page 384
%
%	mmai, 02/08/98
%	--------------

Mo = 10.^(1.5*Mw + 9.05);
