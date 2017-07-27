function [spacing, aperture, density] = coneSize(ecc,ang, varargin)
% [spacing, aperture, density] = coneSize(ecc,ang, varargin)
%
% Calculate expected cone spacing and aperture size at this eccentricity and angle.
%
%  [spacing, aperture, density] = coneSize(ecc,ang, varargin)
%
% Inputs
%  ecc - eccentricity in meters
%  ang - angle in deg (0 on x axis, 90 on y axis)
%  density - cones per mm2
%
% Key/value pairs
%  'whichEye' - 'left','right' [default 'left'] - which eye to compute for.
%     Passed into coneDensity.
%
% Returns
%  spacing - center to center spacing in meters
%  aperture - inner segment linear capture size in meters.  Typically, we
%    set the photoPigment pdHeight and pdWidth both equal to this.
%
% By default, the aperature is set to 0.7*spacing.  We are not sure this is
% a perfect number.
%
% See also: coneDensity.
%
% BW ISETBIO Team, 2016

p = inputParser;
vFunc = @(x)(isnumeric(x) && all(0 <= x & x < 30*1e-3));  % Meters
p.addRequired('ecc',vFunc);
vFunc = @(x)(isnumeric(x) && all(0 <= x <= 3600));    % Angle in degrees
p.addRequired('ang',vFunc);
vFunc = @(x)(ismember(x,{'left','right'}));
p.addParameter('whichEye','left',vFunc);

p.parse(ecc,ang,varargin{:});

ecc = p.Results.ecc;
ang = p.Results.ang;
whichEye = p.Results.whichEye;

% cones/mm2
density = coneDensity(ecc,ang,whichEye);
conesPerMM = sqrt(density);
conesPerM = conesPerMM*1e3;

% Made up for now
spacing = 1./conesPerM;
aperture = 0.7*spacing;  

end