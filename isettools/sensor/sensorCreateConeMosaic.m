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
%   sensor = sensorCreate('human',[],coneP);
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

density = coneGet(coneP,'spatial density');
if isfield(coneP,'rSeed'), rSeed = coneP.rSeed;
else rSeed = []; end

coneAperture  = sensorGet(sensor,'pixel size');   % This should get added to coneCreate/Set/Get
sz      = sensorGet(sensor,'size');

switch ieParamFormat(coneGet(coneP,'species'))
    case 'human'
        % Create a model human sensor array structure.
        sensor = sensorSet(sensor,'name', ...
                        sprintf('human-%.0f',vcCountObjects('ISA')));

        % Deal with the case of black pixels in here
        [xy, coneType] = humanConeMosaic(sz,density,coneAperture(1)*1e6,rSeed);
        coneType = reshape(coneType,sz);
        
        % The cone type defines each cone in CFA, so the size must match
        sensor = sensorSet(sensor,'size',size(coneType));

        % Set the cone pattern (full size)
        sensor = sensorSet(sensor,'pattern',coneType);
        
        % Adjust the spectra so that the first is black and the remaining
        % three are the Stockman fundamentals.
        %
        % The Stockman functions are based on color-matching with the units
        % of the light in energy. ISET computes the sensor response based
        % on an irradiance in photons.  So we use the form of the Stockman
        % filters that is appropriate for photon input.  See the script
        % stockmanQuanta in the data/human directory.
        % fname = fullfile(isetRootPath,'data','human','stockmanQuanta.mat');
        % [fsQuanta,fN] = ieReadColorFilter(wave,fname);
        
        % Instead of loading stockman cone quanta fundamentals, we create a
        % cone structure and compute all ocular transmittance and cone
        % absorptance
        sensor = sensorSet(sensor, 'human cone', coneP);
        fN = {'kBlack', 'rLong', 'gMiddle', 'bShort'};
        sensor = sensorSet(sensor,'filter names',fN);
              
        sensor = sensorSet(sensor,'cone locs',xy);
        sensor = sensorSet(sensor,'cone type',coneType);
        sensor = sensorSet(sensor,'rSeed',rSeed);
               
    otherwise
        error('Unknown species.');
end

return

