function cone = coneCreate(species, varargin)
% Create a cone structure
%   cone = coneCreate(species, varargin)
%
% create a default cone structure, which includes cone absorptance, spatial
% density, pigment density, etc.
%
% Inputs:
%   species   - only support 'human' at this point
%   varargin  - optional input, should be key-value pairs to set some
%               parameters in cone structure
%
% Output:
%   cone      - cone structure
%
% Useful Notes:
%   Absorbance spectra:   Normalized to a peak value of 1.
%   Optical density (OD): Pigment specific
%   Absorptance spectra:  Proportion of quanta actually absorbed
%   The relationship between absorbance and absorptance is as below:
%     absorptance = 1 - 10.^(-OD * absorbance)
%     absorbance  = log10(1 - absorptance) * diag(-1 ./ opticalDensity);
%
% Example:
%    cone = coneCreate('human');
%
% See also:
%   coneGet, coneSet, sensorGet, sensorSet
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check 
if notDefined('species'),     species = 'human'; end
species = ieParamFormat(species);

%% Create cone structure
switch species
    case {'human','humanfovea'}
        
        % Set basic information
        cone.species = 'human';
        cone.type    = 'cone';
        cone.name    = 'default human cone';
        cone.wave    = (400:10:700)';
        
        % photopigment optical densities for L,M,S
        cone.opticalDensity = [0.5 0.5 0.4]';
        
        % load absorbance data
        cone.absorbance = 10.^ieReadSpectra('coneAbsorbance',cone.wave);
        
        % Peak absorptance efficiency
        % Note that the actual peak efficiency is the product of this value
        % and the opticalDensity. The effective value should be around .3
        % for all three types of cones
        cone.peakEfficiency = [2 2 2]/3;

        cone.spatialDensity = [0 .6 .3 .1];  % No blanks, L,M,S
        
    otherwise
        error('Unknown type encountered');
end


%% Over-write default parameters with what user sent in
% Check that we have an even number of arguments
n = length(varargin);
if isodd(n), error('Must have param,val pairs'); end

% Set the parameter,val pairs into the cone structure
for ii=1:2:(n-1)
    cone = coneSet(cone,varargin{ii},varargin{ii+1});
end

end