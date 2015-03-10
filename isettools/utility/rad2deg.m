function d = rad2deg(r,varargin)
% Convert radians to degrees
%
%  d = rad2deg(r)
%
% Also converts to minutes by extra argument.
%   rad2deg(r);
%   rad2deg(r,'arcmin');
%   rad2deg(r,'arcsec');
%
% Copyright ImagEval Consultants, LLC, 2005.

if isempty(varargin), d = (180/pi)*r; 
else                  d = r*ieUnitScaleFactor(varargin{1});
end   

end