function val = coneGet(cone, param, varargin)
% Get cone-specific properties
%
%    val = coneGet(cone, param, varargin)
%
%  The cone properties do not include inert pigments or lens transmittance.
%  Those are accounted for as part of the sensor structure (sensorGet).
%  The information about the lens macular and cone are stored in the slot
%  sensor.human.
%
%  Thus, to determine the spectral quantum efficiency of the cones in the
%  context of the eye, use sensorGet. 
%
%  Inputs:
%    cone     - cone structure, created by coneCreate
%    param    - parameter name to get, spaces and case insensitive
%    varargin - possible units for some parameters
%
%  Outputs:
%    val      - value for parameter, if not found, return empty
%
%  Supported params:
%    {'name'}                        - user defined name of the cone
%    {'type'}                        - should be 'cone'
%    {'species'}                     - cone species, generally 'human'
%    {'cone spatial density'}        - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'wave', 'wavelength'}          - wavelength of samples in cones
%    {'PODs','POD'}                  - pigment density vector for [L,M,S]
%    {'LPOD'}                        - L pigment density
%    {'MPOD'}                        - M pigment density
%    {'SPOD'}                        - S pigment density
%    {'peak efficiency'}             - quantal efficiency
%    {'absorbance'}                  - cone absorbance
%    {'absorptance'}                 - cone absorptance without any
%                                      pre-receprocal transimitance
%    {'quantal fundamentals'}        - quantal fundamentals of the cones
%
%
%  Example:
%    cone = coneCreate('human');
%    expTime = coneGet(cone, 'species');
%
%  See also:
%    coneSet, coneCreate,  sensorGet(sensor,'spectral qe');
%
%  TODO:
%    For most parameters, we should accept a third parameter as wavelength
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('Cone structure required'); end
if notDefined('param'), error('Parameter name required'); end

%% Get property value
param = ieParamFormat(param);  % Lower case and remove spaces
switch param   
    case {'name'}
        val = cone.name;
    case {'type', 'visualfield'}
        val = cone.type; % Should always be 'cone'
    
    case {'wave', 'wavelength'}
        val = cone.wave;
    case {'species', 'kind'}
        % Animal species.
        % Currently supports only 'human'
        val = cone.species;
        
    case {'spatialdensity', 'conespatialdensity'}
        % Spatial density of the cone samples.
        % The vector should be in form [K, L, M, S]
        % Sum of the four elememts should be 1
        val = cone.spatialDensity;
            
    case {'pods','pod'}
        % pigment density inside each cone cell
        % 3-element vector, for [L,M,S] repectively
        val = cone.opticalDensity;
    case {'lpod'}
        val = cone.opticalDensity(1);
    case {'mpod'}
        val = cone.opticalDensity(2);
    case {'spod'}
        val = cone.opticalDensity(3);
        
    case {'peakefficiency'}
        val = cone.peakEfficiency;
        
    case {'absorbance'}
        val = cone.absorbance;
        
    case {'conespectralabsorptance','absorptance'}
        % coneGet(cone,'cone spectral absorptance')
        %
        % This is the cone absorptance without the ocular media
        absorbance = coneGet(cone, 'absorbance');
        PODs = coneGet(cone, 'PODs');
        val = 1-10.^(-absorbance*diag(PODs));
        
    case {'quantalfundamentals','photonfundamentals'}
        % coneGet(cone,'photon fundamentals')
        %
        % Cone absorptance scaled by peak efficiency
        % Note that this is for pure cones, without any oclus absorptance
        val = coneGet(cone, 'cone spectral absorptance');
        qe  = coneGet(cone, 'peak efficiency');
        if length(qe) == size(val,2)
            for ii = 1 : size(val, 2)
                val(:,ii) = val(:,ii) * qe(ii);
            end
        end
        val = val ./ repmat(max(val), [size(val, 1) 1]);
        
    case {'energyfundamentals'}
        % coneGet(cone, 'energy fundamentals')
        %
        % cone absorptance in energy
        % Note that this is for pure cones, without any oclus absorptance
        val = coneGet(cone, 'quantal fundamentals');
        % Get constants
        h = vcConstants('planck');
        c = vcConstants('speed of light');
        wave = coneGet(cone, 'wave');
        % Convert to energy. Note that the sensitivity is the inverse of
        % spetral absorption. Thus, we use conversion of energy to quanta
        % on quantal fundamentals
        val = val / (h*c) .* (1e-9 * repmat(wave, [1 size(val, 2)])); 
        
        % Normalize
        val = val ./ repmat(max(val), [size(val, 1) 1]);
        
    otherwise
        error('Unknown parameter encountered');
end

end
