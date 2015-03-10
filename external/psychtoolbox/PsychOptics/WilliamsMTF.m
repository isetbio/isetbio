function mtf = WilliamsMTF(s)
% mtf = WilliamsMTF(s)
% 
% Compute the MTF measured by Williams et. al.
% for 633 nm light.
%
% 7/11/94		dhb		Wrote it.
% 7/14/94		dhb		Pulled calculation of rest into separate function.

diff = DiffractionMTF(s,3,633);
rest = WilliamsRestMTF(s);
mtf = diff .* rest;
