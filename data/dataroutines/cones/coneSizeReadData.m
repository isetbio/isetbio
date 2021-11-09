function [spacing, aperture, density, params, comment] = coneSizeReadData(varargin)
%%coneSizeReadData  Read in data about cone size parameters
%
% Syntax:
%    [spacing, aperture, density] = coneSizeReadData;
%
% Descirption:
%     Calculate expected cone spacing and aperture size at this eccentricity and angle.
%     This is done based on cone density, obtained via parameter coneDensitySource.
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
% Optional key/value pairs
%    'species'                  What species?
%                                 'human' (default)
%
%    'coneDensitySource'        Source for cone density estimate, on which other values are based.
%                               This is passed on to coneDensityReadData.  See help for that function.
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
% See also: coneDensityReadData.

% BW ISETBIO Team, 2016
%
% 08/16/17  dhb  Call through new coneDensityReadData rather than old coneDensity.
% 02/17/19  npc  Added useParfor k/v pair

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
p.addParameter('useParfor', false, @(x)(((isempty(x))||islogical(x))));
p.parse(varargin{:});

%% Set up params return.
params = p.Results;

%% Take care of case where a function handle is specified as source
%
% This allows for custom data to be defined by a user, via a function that
% could live outside of ISETBio.
%
% This function needs to handle 
if (isa(params.coneDensitySource,'function_handle'))
    [spacing, aperture, density, comment] = params.coneSizeSource(varargin{:});
    return;
end

%% Get density.  This can just take the params structure, except we change the source name
[density,~,comment] = coneDensityReadData(varargin{:});
conesPerMM = sqrt(density);
conesPerM = conesPerMM*1e3;

%% Compute spacing and aperture
spacing = 1./conesPerM;
aperture = 0.7*spacing;  

end