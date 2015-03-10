function cone = coneSet(cone, param, val, varargin)
% Set parameters in a cone structure.
%
%      cone = coneSet(cone, param, val, [varargin])
%
%  Inputs:
%    cone     - cone structure, created by coneCreate
%    param    - parameter name to set, spaces and case insensitive
%    val      - new value of parameter to be set to
%    varargin - possible units for some parameters
%
%  Outputs:
%    cone     - cone structure with parameter and its related information
%               set to new value
%
%  Supported parameters:
%    {'name'}                        - name of the cone structure
%    {'spatial density'}             - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'wave', 'wavelength'}          - wavelength of samples in cones
%    {'lens'}                        - eye lens structure
%    {'lens density'}                - lens OD
%    {'macular'}                     - macular pigment structure
%    {'macular density'}             - macular density
%    {'PODs','POD'}                  - PODs vector for [L,M,S]
%    {'LPOD'}                        - L POD density
%    {'MPOD'}                        - M POD density
%    {'SPOD'}                        - S POD density
%    {'peak lambda', 'lambda max'}   - peak spectra position, nomogram
%    {'qe', 'peak efficiency'}       - quantal efficiency
%    {'absorbance'}                  - cone absorbance, not recommended to
%                                      set this parameter directly unless
%                                      you know what you are doing
%    {'adaptation type','adapt type'}- cone adaptation type, see
%                                      coneAdaptation for more details
%
%    MORE SUPPORTED PARAMETERS CAN BE FOUND IN FUNCTION sensorSet
%
%  Example:
%    cone = coneCreate('human');
%    cone = coneSet('spatial density',[0.1 0.65 0.2 0.05]);
%
%  See also:
%    coneCreate, coneGet, sensorSet
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('cone structure required'); end
if notDefined('param'), error('parameter name required'); end
if notDefined('val'), error('new value of parameter required'); end

%% Set parameters
param = ieParamFormat(param);  % Lower case and remove spaces
switch param
    case {'name'}
        cone.name = val;
    case {'spatialdensity', 'conedensity'}
        val = val(:);
        if length(val) == 3 % Cone density given in [L,M,S] format
            val = [1-sum(val); val];
        else
            assert(length(val)==4, 'Unknown density format encountered');
        end
        cone.spatialDensity = val;
    
    case {'wave', 'wavelength'}
        % Make sure wavelengths are consistent.
        % We interpolate the unit density data in lens and macular by
        % calling lensSet(lens, 'wave') and macularSet(...)
        val = val(:);
        cone.wave = val;
        cone.macular = macularSet(cone.macular, 'wave', val);
        cone.lens    = lensSet(cone.lens, 'wave', val);

    case {'lens'}
        % Set the lens structure
        % We should first make sure that val is actually a lens structure
        % and it has the same sampling wavelength as the cone structure
        assert(strcmp(val.type, 'lens'), 'val should be lens structure');
        cone.lens = lensSet(val, 'wave', coneGet(cone, 'wave'));
        
    case {'lensdensity'}
        % lens OD
        cone.lens = lensSet(cone.lens, 'density', val);
        
    case {'macular'}
        % Set macular structure
        % We should first make sure that val is actually a macular
        % structure and it has the same wavelength as the cone
        assert(strcmp(val.type, 'macular'),'val should be macular struct');
        cone.macular = macularSet(val, 'wave', coneGet(cone, 'wave'));
        
    case {'maculardensity'}
        % cone = coneSet(cone,'macular density',val)
        % val is typically between 0 and 0.7, a range of macular pigment
        % densities.
        m    = coneGet(cone,'macular');
        m    = macularSet(m,'density',val);
        cone = coneSet(cone,'macular',m);
        
    case {'pods','pod'}
        % Pigment optical densities for the cones
        if (any(size(cone.PODs)~= size(val)))
            error('PODs size mismatch');
        end
        cone.PODs = val;
        
    case {'lpod'}
        if ~isscalar(val), error('val should be scalar'); end
        cone.PODS(1) = val;
    case {'mpod'}
        if ~isscalar(val), error('val should be scalar'); end
        cone.PODS(2) = val;
    case {'spod'}
        if ~isscalar(val), error('val should be scalar'); end
        cone.PODS(3) = val;
        
        %     case {'peaklambda', 'lambdamax'}
        %         cone.peakShift = val;
        %         cone.absorbance = StockmanSharpeNomogram(cone.wave, val)';
        %         cone.absorbance = padarray(cone.absorbance,[0 1],'pre');
        %
    case {'peakefficiency', 'qe', 'quantalefficiency'}
        % Peak absorptance
        cone.peakEfficiency = val;
        
    case {'absorbance'}
        disp('Overwrite absorbance, please be sure what you are doing');
        cone.absorbance = val;
    
    case {'adaptationtype', 'adapttype'}
        if cone.adaptType ~= val
            cone = coneAdapt(cone, val);
        end
        
    otherwise 
        error('Unknown parameter encountered');
end

end
