function dat = replaceNaN(dat,val)
% Replace NaNs with val
%
%  dat = replaceNaN(dat,val)
%
% Copyright ImagEval Consultants, LLC, 2005.

warning('Deprecated. Use dat(isnan(dat))= val; instead');
l = isnan(dat);
dat(l) = val;

end