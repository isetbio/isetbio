function oiSum = oiAdd(oi1,oi2,wgts)
% We need to write a real and tested function here.
% But for now ... we need this as a quick and dirty test

oiSum = oi1;
oiSum.data.photons = wgts(1)* oi1.data.photons + wgts(2)*oi2.data.photons;

end