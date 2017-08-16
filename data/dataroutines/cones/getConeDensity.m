function [coneDensity,params,comment] = getConeDensity(varargin)
%%getConeDensity  Compute cone packing density as a function of retinal position
%
% Syntax:
%     [coneDensity,params,comment] = getConeDensity;
%     coneDensity = getConeDensity('eccentricity',8*1e-3,'angle',10,'whichEye','left');
%
% Description:
%     Compute cone packing density as a function of retinal position.
%
%     For the left eye, the coordinate system is 0 degrees angle corresponds
%     to the nasal meridian, 90 degrees to superior, 180 to temporal and 270
%     to inferior.
%
%     For the right eye, the coordinate system is 0 degrees angle corresponds
%     to the temporal meridian, 90 degrees to superior, 180 to nasal, and 270
%     to inferior.
%
% Input:
%     None.
%
% Output:
%     coneDensity                Cone packing density in cones/mm^2
%
%     params                     Structure of key/value pairs used to generate data.
%
%     comment                    A short comment describing the data, returned as a string.
%
% Optional key/value pairs
%    'species'                  What species?
%                                 'human' (default)
%
%    'coneDensitySource'        Source for cone density estimate
%                                 'Curcio1990' (default). From Figure 6 of Ref 1 below.  This is in data file coneDensityCurcio1990Fig1.
%
%                                  The value for 'coneDensitySource' may be passed as a function handle, in
%                                  which case the passed function is called direclty with the key/value pairs passed to this
%                                  routine. The passed function must return the same values as getConeDensity does.
%
%    'eccentricity'             Retinal eccentricity in meters of retina, default is 0
%                               There are 1000 mm/m and about 0.30 mm/deg.
%
%    'angle'                    Polar angle of retinal position in degrees (default 0).
%
%    'whichEye'                 Which eye, 'left' or 'right' (default 'left').
%
% References:
%   1) Curcio, C. A., Sloan, K. R., Kalina, R. E. and Hendrickson, A. E.
%      (1990), Human photoreceptor topography. J. Comp. Neurol., 292:
%      497?523. doi: 10.1002/cne.902920402
%   2) Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., & Burns, S. A.
%      (2011). Variation of Cone Photoreceptor Packing Density with Retinal
%      Eccentricity and Age. Investigative Ophthalmology & Visual Science,
%      52(10), 7376?7384. http://doi.org/10.1167/iovs.11-7199
%
% See also: coneMosaic, coneSize.

% HJ, ISETBIO TEAM, 2015

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species','human', @ischar);
p.addParameter('coneDensitySource','Curcio1990',@(x) (ischar(x) | isa(x,'function_handle')));
p.addParameter('eccentricity',0, @isnumeric);
p.addParameter('angle',0, @isnumeric);
p.addParameter('whichEye','left',@ischar);
p.parse(varargin{:});

%% Set up params return.
params = p.Results;

%% Take care of case where a function handle is specified as source
%
% This allows for custom data to be defined by a user, via a function that
% could live outside of ISETBio.
if (isa(params.coneDensitySource,'function_handle'))
    [coneDensity,params,comment] = params.coneDensitySource(varargin{:});
    return;
end

%% Handle choices
switch (params.species)
    case {'human'}
        switch (params.coneDensitySource)
            case 'Curcio1990'
                % Load the digitized cone density from the ISETBio style mat file.  The
                % data file has separate structs for inferior, nasal, superior and temporal meridians.
                % These each have fields 'density' as a function of 'eccMM' in units of cones/mm2.
                theData = getRawData('coneDensityCurcio1990Fig6','datatype','isetbiomatfileonpath');
                
                % Convert eccentricity from meters to mm
                eccMM = params.eccentricity*0.3;
                
                % Set up for interpolation for retinal position amplitude on each axis (nasal, superior,
                % temporal and inferior)
                onAxisD = zeros(5, numel(params.eccentricity));
                angleQ = [0 90 180 270 360];
                
                % Compute packing density for superior and inferior
                onAxisD(2,:) = interp1(theData.superior.eccMM, theData.superior.density, eccMM);
                onAxisD(4,:) = interp1(theData.inferior.eccMM, theData.inferior.density, eccMM);
                
                % Nasal and temporal, respecting our le/re coordinate convention.
                switch lower(params.whichEye)
                    case 'left'
                        onAxisD(1,:) = interp1(theData.nasal.eccMM, theData.nasal.density, eccMM);
                        onAxisD(3,:) = interp1(theData.temporal.eccMM, theData.temporal.density, eccMM);
                    case 'right'
                        onAxisD(1,:) = interp1(theData.temporal.eccMM, theData.temporal.density, eccMM);
                        onAxisD(3,:) = interp1(theData.nasal.eccMM, theData.nasal.density, eccMM);
                    otherwise
                        error('unknown input for whichEye');
                end
                onAxisD(5,:) = onAxisD(1,:);
                
                % Interpolate for angle
                coneDensity = interp1(angleQ, onAxisD, params.angle, 'linear');
                
                comment = 'Cone density derived from Figure 6 of Curcio et al (1990).  See getConeDensity.';
                
            otherwise
                error('Unsupprted source specified');
        end
        
    otherwise
        error('Unsupported species specified');
end

end
