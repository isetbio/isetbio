function oiSum = oiAdd(oi1,oi2,wgts)
% OIADD - Add the photons in two matched oi structs
%
% We use this to form the spectral irradiance images in an @oiSequence. 
%
%   oiSum = oiAdd(oi1,oi2,wgts)
%
% We should check that the two oi structs match, and then we add the
% photons from the second oi into the first oi.
%
% Example
%    wgts = [0.5, 0.5]; 
%    oiSum = oiAdd(oi1,oi2,wgts);
%
% See also OISEQUENCE
%
% BW ISETBIO Team, 2016

oiSum = oi1;

oiSum.data.photons = wgts(1)* oi1.data.photons + wgts(2)*oi2.data.photons;

end