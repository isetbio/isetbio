function [spacing, aperture, density] = coneSizeReadData(ecc,ang, varargin)
%%coneSizeReadData  Read in data about cone size parameters
%
% Syntax:
%    [spacing, aperture, density] = coneSizeReadData(ecc,ang)
%
% Descirption:
%     Calculate expected cone spacing and aperture size at this eccentricity and angle.
%
%     The coordinate system is as defined by coneDensityReadData.
%
%     By default, the aperature is set to 0.7*spacing.  We are not sure this is
%     a perfect number.
%
% Input:
%     ecc          Eccentricity in meters.
%
%     ang          Angle in degrees.
%
% Output:
%     spacing      Center to center spacing in meters.
%
%     aperture     Inner segment linear capture size in meters.  Typically, we
%                  set the photoPigment pdHeight and pdWidth both equal to this.
%
%     density      Cones per mm2. This is the density returned by coneDensityReadData.
%
% Optional key/value pairs:
%  'whichEye' - 'left','right' [default 'left'] - which eye to compute for.
%     Passed into coneDensity.
%
% See also: coneDensityReadData.

% BW ISETBIO Team, 2016
%
% 08/16/17  dhb  Call through new coneDensityReadData rather than old coneDensity.

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
density = coneDensityReadData('eccentricity',ecc,'angle',ang,'whichEye',whichEye);
conesPerMM = sqrt(density);
conesPerM = conesPerMM*1e3;

% Made up for now
spacing = 1./conesPerM;
aperture = 0.7*spacing;  

end