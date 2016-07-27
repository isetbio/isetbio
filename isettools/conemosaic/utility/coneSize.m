function [spacing, aperture, density] = coneSize(ecc,ang, varargin)
% Calculate expected cone spacing and aperture size at this eccentricity.
%
%  [spacing, aperture, density] = coneSize(ecc,ang, varargin)
%
% Inputs
%  ecc - eccentricity in meters
%  ang - angle in deg (0 on x axis, 90 on y axis)
%' density - cones per mm2
%
% Returns
%  spacing - center to center spacing in meters
%  aperture - inner segment capture size in meters
%
% Dummy routine to get cone spacing and size at different eccentricities.
% We will set the photoPigment pdSize,pdWidth based on this.  We will set
% the cone size and width
%
% BW ISETBIO Team, 2016

p = inputParser;
vFunc = @(x)(isnumeric(x) && 0 <= x && x < 30*1e-3);  % Meters
p.addRequired('ecc',vFunc);
vFunc = @(x)(isnumeric(x) && 0 <= x <= 360);    % Angle in degrees
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
spacing = 1/conesPerM;
aperture = 0.7*spacing;   % Rods .... need to get the right value for this

end