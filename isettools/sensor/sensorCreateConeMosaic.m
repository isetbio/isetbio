function sensor = sensorCreateConeMosaic(sensor,coneP)
%Create a sensor with a random cone mosaic
%
%  sensor = sensorCreateConeMosaic(sensor, coneP)
%
% The human mosaic is a random array of empty (K), L, M, and S cones. THe
% parameters of the coneP (cone structure) object defines many features of
% the cone mosaic, including the relative cone densities, peak
% efficiencies, absorbance, and optical Density.
%
%
% Inputs
%  sensor: Initialized sensor
%  params: cone structure (see coneCreate)
%    .sz:            72,88 mosaic (default)
%    .density:     
%       human:  [0 0.6, 0.3, 0.1] default
%       mouse:  Not yet implemented
%           Any two non-zero numbers work for the mouse. A zero for one
%           cone type will create a monochromatic sensor of the other cone
%           type.
%    .coneAperture:  [1.5 1.5]*1e-6 (default). Microns. Probably too small.
%    .rSeed:         Random number seed for creating the mosaic
%
% Returns
%  sensor:   Human sensor
%  xy:       Cone xy positions
%  coneType: Vector of cone types 1:4 for K,L,M,S
%
% See also:  For an alternative method use:
%
%   coneP = coneCreate;  % Then adjust the coneP values using coneSet/Get
%   sensor = sensorCreate('human', coneP);
%
% Examples:
%  sensor = sensorCreate('human');
%  sensor = sensorCreateConeMosaic(sensor,coneCreate);
%  xy = sensorGet(sensor, 'cone xy'); 
%  coneType = sensorGet(sensor, 'coneType');
%
%  vcNewGraphWin; plot(sensorGet(sensor,'wave'),sensorGet(sensor,'spectralQE'));
%  vcNewGraphWin; conePlot(xy,coneType);
%  vcNewGraphWin; sensorConePlot(sensor)
%
% See also:  sensorCreate('human'), coneCreate
%
% (c) Copyright, 2010, ImagEval

if notDefined('sensor'), error('Human sensor required'); end

d = coneGet(coneP,'spatial density');
if isfield(coneP, 'rSeed'), rSeed = coneP.rSeed;
else rSeed = []; end

coneWidth = sensorGet(sensor, 'pixel deltax', 'um');
sz = sensorGet(sensor, 'size');

switch ieParamFormat(coneGet(coneP,'species'))
    case 'human'
        % Create a model human sensor array structure
        name = sprintf('human-%.0f',vcCountObjects('sensor'));
        sensor = sensorSet(sensor,'name', name);

        % Deal with the case of black pixels in here
        [xy, coneType] = humanConeMosaic(sz, d, coneWidth, rSeed);
        coneType = reshape(coneType, sz);
        
        % Apply the cone structure to sensor
        sensor = sensorSet(sensor, 'human cone', coneP);
        
        % Set filter names
        fN = {'kBlack', 'rLong', 'gMiddle', 'bShort'};
        sensor = sensorSet(sensor, 'filter names', fN);
        
        % Set cone type and relative positions
        sensor = sensorSet(sensor,'cone locs',xy);
        sensor = sensorSet(sensor,'cone type',coneType);
        sensor = sensorSet(sensor,'rSeed',rSeed);            
    otherwise
        error('Unknown species.');
end

end