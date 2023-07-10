% t_conesCustomMosaic
%
% Show how to customize various mosaic parameters
%
% Description:
%    Demonstrates how to customize various parameters associated with the
%    cone mosaic.
%

%% Jon W.'s desired coordinates
eccentricityDeg = 6;
angle = 90;

%% Convert to x, y position in meters
eccentricityM = 1e-3 * 0.3 * eccentricityDeg;
eccentricityXM = cosd(angle) * eccentricityM;
eccentricityYM = sind(angle) * eccentricityM;

%% Create coneMosaic using 'Song2011Young' density data
% Pigment width/height and pdWidth/pdHeight should vary with data source
cm = coneMosaic('center', [eccentricityXM eccentricityYM], ...
    'coneDensitySource', 'Song2011Young');
cm.pigment

%% So instead try the old subject data to find out
cm = coneMosaic('center', [eccentricityXM eccentricityYM], ...
    'coneDensitySource', 'Song2011Old');
cm.pigment

%% Make sure passing units parameter doesn't screw things up for now.

% Eventually we should respect all of the parameters handled by
% coneDensityReadData, but for right now center must be passed in meters.
% We'll fix this when we update cone mosaic to use eccentricity and angle.
cm = coneMosaic('center', [eccentricityXM eccentricityYM], ...
    'coneDensitySource', 'Song2011Old', 'eccentricityUnits', 'mm');
cm.pigment

%%