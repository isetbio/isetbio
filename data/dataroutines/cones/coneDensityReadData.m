function [coneDensity,params,comment] = coneDensityReadData(varargin)
% Compute cone packing density as a function of retinal position
%
% Syntax:
%     [coneDensity,params,comment] = coneDensityReadData;
%     coneDensity = coneDensityReadData('eccentricity',8*1e-3,'angle',10,'whichEye','left');
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
%     If an eccentricity beyond the range of the underlying data source is requested,
%     NaN is returned.
%
% Input:
%     None.
%
% Output:
%     coneDensity                Cone packing density in cones/mm^2. Returned as a row vector
%                                of same length as eccentricity and angle.
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
%                                 'Curcio1990'         From Figure 6 of Ref 1 below (default).
%                                 'Song2011Old'        From Table 1 of Ref 2 below, old subjects data.
%                                 'Song2011Young'      From Table 1 of Ref 2 below, young subjects data.
%
%                                  The value for 'coneDensitySource' may be passed as a function handle, in
%                                  which case the passed function is called direclty with the key/value pairs passed to this
%                                  routine. The passed function must return the same values as coneDensityReadData does, and handle
%                                  the same key/value pairs (including unit specification).
%
%    'eccentricity'             Retinal eccentricity, default is 0.  Units according
%                               to eccentricityUnits.  May be a vector,
%                               must have same length as angle.
%
%    'angle'                    Polar angle of retinal position in degrees (default 0).  Units
%                               according to angleUnits.  May be a vector, must have
%                               same size as eccentricity.
%
%    'whichEye'                 Which eye, 'left' or 'right' (default 'left').
%
%    'eccentriticyUnits'        String specifying units for eccentricity.
%                                  'm'                  Meters (default).             
%                                  'mm'                 Millimeters.
%                                  'um'                 Micrometers.
%                                  'deg'                Degrees of visual angle, 0.3 mm/deg.
%
%    'angleUnits'               String specifying units  for angle.
%                                  'deg'                Degrees (default).
%                                  'rad'                Radians.
%
%    'useParfor'                Logical. Default false. Used to parallelize
%                               the interp1 function calls which take a
%                               long time. This is useful when generating
%                               large > 5 deg mosaics.
%
% References:
%   1) Curcio, C. A., Sloan, K. R., Kalina, R. E. and Hendrickson, A. E.
%      (1990), Human photoreceptor topography. J. Comp. Neurol., 292:
%      497?523. doi: 10.1002/cne.902920402
%   2) Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., & Burns, S. A.
%      (2011). Variation of Cone Photoreceptor Packing Density with Retinal
%      Eccentricity and Age. Investigative Ophthalmology & Visual Science,
%      52(10), 7376-7384. http://doi.org/10.1167/iovs.11-7199
%
% See also: coneMosaic, coneSize, coneDensitySong2011MakeData.

% HJ, ISETBIO TEAM, 2015
%
% 08/16/17  dhb  Big rewrite.
% 02/17/19  npc  Added useParfor k/v pair
%
%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species','human', @ischar);
p.addParameter('coneDensitySource','Curcio1990',@(x) (ischar(x) | isa(x,'function_handle')));
p.addParameter('eccentricity',0, @isnumeric);
p.addParameter('angle',0, @isnumeric);
p.addParameter('whichEye','left',@(x)(ismember(x,{'left','right'})));
p.addParameter('eccentricityUnits','m',@ischar);
p.addParameter('angleUnits','deg',@ischar);
p.addParameter('useParfor', false, @(x)((islogical(x))||(isempty(x))));
p.parse(varargin{:});

%% Set up params return.
params = p.Results;

if (isempty(params.useParfor))
    params.useParfor = false;
end

%% Take care of case where a function handle is specified as source
%
% This allows for custom data to be defined by a user, via a function that
% could live outside of ISETBio.
%
% This function needs to handle 
if (isa(params.coneDensitySource,'function_handle'))
    [coneDensity,params,comment] = params.coneDensitySource(varargin{:});
    return;
end


%% Handle units
switch (params.eccentricityUnits)
    case 'm'
        eccMM = params.eccentricity*1e3;
    case 'mm'
        eccMM = params.eccentricity;
    case 'um'
        eccMM = params.eccentricity*1e-3;
    case 'deg'
        eccMM = params.eccentricity*0.3;
    otherwise
        error('Unknonwn units for eccentricity specified');
end
params.eccentricity = NaN;

switch (params.angleUnits)
    case 'deg'
        angleDeg = params.angle;
    case 'rad'
        angleDeg = 180*params.angle/pi;
    otherwise
        error('Unkonwn units for angle specified');
end
params.angle = NaN;

%% Make sure eccentricity and angle are vectors of the same length.
if (~isvector(eccMM) || ~isvector(angleDeg))
    error('Passed eccentricity and angle arguments must be vectors.');
end
if (length(eccMM) ~= length(angleDeg))
    error('Passed eccentricity and angle vectors must have the same length.')
end

%% Wrap angles so that they are in range 0 <= angle < 360
index = find(angleDeg < 0);
while (~isempty(index))
    angleDeg(index) = angleDeg(index)+360;
    index = find(angleDeg < 0);
end
index = find(angleDeg >= 360);
while (~isempty(index))
    angleDeg(index) = angleDeg(index)-360;
    index = find(angleDeg >= 360);
end

%% Handle choices
switch (params.species)
    case {'human'}
        switch (params.coneDensitySource)
            case {'Curcio1990', 'Song2011Old', 'Song2011Young'}
                % Load the digitized cone density from the ISETBio style mat file.  The
                % data file has separate structs for inferior, nasal, superior and temporal meridians.
                % These each have fields 'density' as a function of 'eccMM' in units of cones/mm2.
                switch (params.coneDensitySource)
                    case 'Curcio1990'
                        theData = rawDataReadData('coneDensityCurcio1990','datatype','isetbiomatfileonpath');
                    case 'Song2011Old'
                        theData = rawDataReadData('coneDensitySong2011Old','datatype','isetbiomatfileonpath');
                    case 'Song2011Young'
                        theData = rawDataReadData('coneDensitySong2011Young','datatype','isetbiomatfileonpath');
                end
                
                % Set up for interpolation for retinal position amplitude on each axis (nasal, superior,
                % temporal and inferior)
                onAxisD = zeros(5, numel(eccMM));
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
                
                % Interpolate for each angle
                coneDensity = zeros(1,length(eccMM));
                
                if (params.useParfor)
                    % For large mosaics (>5 deg), having a parfor loop here
                    % reduced computation time significantly because this
                    % function is called all the time during the iterative
                    % process.
                    parfor aa = 1:length(eccMM)
                        coneDensity(aa) = interp1(angleQ, onAxisD(:,aa), angleDeg(aa), 'linear');
                    end
                else
                    for aa = 1:length(eccMM)
                        coneDensity(aa) = interp1(angleQ, onAxisD(:,aa), angleDeg(aa), 'linear');
                    end
                end
                comment = 'Cone density derived from Figure 6 of Curcio et al (1990).  See coneDensityReadData.';
                
            otherwise
                error('Unsupported source specified');
        end
        
    otherwise
        error('Unsupported species specified');
end

end
