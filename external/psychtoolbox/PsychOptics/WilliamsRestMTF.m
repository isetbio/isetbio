function mtf = WilliamsRestMTF(s)
% mtf = WilliamsRestMTF(s)
%
% Compute the portion of the MTF measured by
% Williams et. al. that is not accounted for
% by diffraction.
%
% 7/14/94		dhb		Wrote it.

a = 0.1212;
w1 = 0.3481;
w2 = 0.6519;
mtf = w1*ones(size(s)) + w2*exp(-a*s);
