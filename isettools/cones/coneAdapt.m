function [sensor, adaptedData] = coneAdapt(sensor, typeAdapt, params)
%% Cone adaptation
%
%   [sensor, adaptedData] = coneAdaptation(sensor, typeAdapt, varargin)
%
% Implement adaptation models to produce the cone volts. Cone absorption
% sampels should be computed and stored in sensorGet(sensor, volts).
%
% The mean level is determined as a reference so that increments and
% decrements are meaningful. This is needed so that on- and off-cells can
% work. When we have a cone absorption that is above the mean, this should
% cause an off-cell to reduce its firing rate and on-cell to increase.
% Conversely, cone absorptions below the mean cause an on-cell to decrease
% and an off-cell to increase.
%
% The discussion in Rodieck book - The first steps in seeing - p. 177 and
% thereabouts for ideas about cone signals. The physiological question is
% what is the meaning of the offset or zero mean. Rodieck describes the
% effect of light as setting the mean transmitter release. In the dark,
% there is a relatively large dynamic range. As the light is on steadily,
% the release decreases and the range available for another flash. If you
% darken from the mean, the rate can increase.
%
% In addition to Rodieck's discussion, there are famous papers by Boynton
% (e.g. adaptation in cones) expressing such a model based on cone ERPs.
%
% The other issue is the total gain.  The cones can only modulate over a
% dynamic range around the mean.  So, we set the gain as well to keep the
% modulation within some range, such as +/- 40 mV of the mean (deplorizing
% voltage).
%
% The third issue is whether all the cones are the same, or there is
% space-variant adaptation. We now here consider the difference between
% different cone types, but we think it's spatial invariant.
%
% Inputs:
%  sensor     - ISETBio sensor with cone absorption computed and stored
%  typeAdapt  - The adaptation model, meaning of value are as below
%           0 = no adaptation
%           1 = a single gain map for the whole cone array
%           2 = one gain map computed for each cone class
%           3 = non-linear adaptation
%           4 = cone adaptation by physiological differential equation
%               (Rieke)
%  params     - Cone adaptation parameters, could include
%    .bgVolts = background voltage
%    .vSwing  = maximum volts response
%
%
% Output:
%   sensor      - ISETBio sensor with cone adaptation parameters set. The
%                 parameters include gain and offset. You could retrieve
%                 the adaptation parameters and data by sensorGet function
%                 calls.
%   adaptedData - Adapted voltage 3D matrix, would be in same size as volts
%                 image in sensor.
%
% Notes:
%   In the future, we will look at the time series and have a time-varying
%   adaptation function. Generally, cones take 200ms and rods take 800ms
%   for adaptation.
%
% Examples:
%   scene = sceneCreate; oi = oiCreate; oi = oiCompute(oi,scene);
%   sensor = sensorCreate('human'); sensor = sensorCompute(sensor,oi);
%
%   sensor = coneAdapt(sensor, 1);
%   adaptedData = sensorGet(sensor, 'adapted volts');
%   gain = sensorGet(sensor, 'adaptation gain');
%   offset = sensorGet(sensor, 'adaptation volts');
%
%   sensor = coneAdapt(sensor,4);  % Rieke model 
%   v = sensorGet(sensor,'volts');
%
% (c) Stanford VISTA Lab, 2014

%% Check inputs and Init
if notDefined('sensor'),      error('sensor is required'); end
if ~sensorCheckHuman(sensor), error('unsupported species'); end
if notDefined('typeAdapt'),   typeAdapt = 3; end
if notDefined('params'),      params = []; end

%% Compute cone adaptations
volts  = double(sensorGet(sensor, 'volts'));
if isempty(volts), error('cone absorptions should be pre-computed'); end

if isfield(params, 'vSwing'),    vSwing = params.vSwing;
else                             vSwing = sensorGet(sensor,'pixel voltageSwing');
end

if isfield(params, 'bgVolts'),     bgVolts = params.bgVolts;
else                               bgVolts = mean(volts(:));
end

if ischar(typeAdapt), typeAdapt = ieParamFormat(typeAdapt); end
    
switch typeAdapt
    case {0, 'noadapt', 'noadaptation'} % no adaptation
        gainMap = 1;
        bgVolts  = 0;
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);
        
    case {1, 'globalgain'}
        % Use same gain for all cone type
        
        % Adjust for bg volts as offset
        volts = volts - bgVolts;
        
        % Set the gain so that the max - min get scaled to vSwig mV
        gainMap = vSwing / (max(volts(:)) - min(volts(:)));
        adaptedData = gainMap * volts;
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);
        
    case {2, 'adaptbyconetype'}
        % Use different gains for each cone type
        % For human, the cone types in sensor are 1~4 for K, L, M, S
        if isscalar(bgVolts), bgVolts = ones(4,1)*bgVolts; end
        assert(numel(bgVolts) == 4, 'bgVolts should be of 4 elements');
        gainMap = ones(4, 1);
        
        % Adjust for backgroud
        volts = volts - bgVolts(sensorGet(sensor, 'cone type'));
        
        % Compute gain map
        for ii = 2 : 4 % L,M,S and we don't need to compute for K
            v = sensorGet(sensor,'volts',ii);
            if ~isempty(v)
                gainMap(ii) = vSwing / (max(v) - min(v));
            end
        end
        
        nSamples = size(volts, 3);
        gainMap = gainMap(sensorGet(sensor, 'cone type'));
        
        adaptedData = volts .* repmat(gainMap, [1 1 nSamples]);
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);

    case {3, 'felice'}
        % In this case, we will do non-linear cone adaptation
        % Reference:
        %   Felice A. Dunn, et al. Light adaptation in cone vision involves
        %   switching between receptor and post-receptor sites,
        %   doi:10.1038/nature06150
        % From there and some other papers, we know two facts:
        %   1) The steady state response follows:
        %      R/R_max = 1/(1 + (I_{1/2}/I_B)^n)
        %   2) The cone adaptation follows
        %      AMP/AMP_{dark} = ((a + bI_B)/(a + I_B))/(cI_B+1)
        % To compute the adapted voltage for input I on background I_B, we
        % use the following process
        %   1) Compute steady state response for I as R_i
        %   2) Compute dark adapted response for I as R_i * AMP / AMP{dark}
        %   3) Compute I_B adapted response
        %
        % Notes:
        %   In the original paper, the formulas are only tested and fitted
        %   for L cones. Here, we use it for all cone types and this might
        %   lead to simulation errors.
        %
        %   Also, by this method, the adapted cone volts are still positive
        %   This is different from the previous two method
        %
        % (HJ) ISETBIO TEAM, 2014
        
        % Check inputs
        if ~isscalar(bgVolts), error('bgVolts should be scalar'); end
        
        % Init function handles
        cg = 1/pixelGet(sensorGet(sensor, 'pixel'), 'conversion gain');
        cg = cg / sensorGet(sensor, 'exp time');
        steadyR  = @(x) 1./(1+(45000 ./ x / cg).^0.7)*vSwing;
        ampRatio = @(x) (100 + 1.3*x/cg)./(100 + x/cg)./(0.00029*x/cg + 1);
        
        % Compute steady state response
        sR = steadyR(volts);
        
        % Compute dark adapted response
        dR = sR ./ ampRatio(volts);
        
        % Compute adapted response at level I_B
        adaptedData = dR * ampRatio(bgVolts);
        
        % Set gain map and offset
        gainMap = adaptedData ./ volts;
        bgVolts = 0; % We don't have an actual offset
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);

    case {4, 'rieke'}
        % See rieke<TAB> for explanations.
        
        % sensor = sensorCreate('human');
        % sensor = sensorSet(sensor,'pixel voltage swing',0.05);
        % v = rand(32,32,200)*sensorGet(sensor,'pixel voltage swing');
        % sensor = sensorSet(sensor,'volts',v);
        % vcAddObject(sensor); sensorWindow;
        % [~,adaptedData] = coneAdapt(sensor, 4);
        % s = 255 / max(abs(adaptedData(:)));
        %                 vcNewGraphWin
        %                 for ii=1:size(adaptedData,3)
        %                     imagesc(abs(s*adaptedData(:,:,ii)));
        %                     pause(0.02);
        %                 end
        % 
        % Fix this mplay(adaptedData,'I');
        
        p = riekeInit;
        expTime = sensorGet(sensor,'exposure time');
        sz = sensorGet(sensor,'size');
        
        % absRate = sensorGet(sensor,'absorptions per second');        
        pRate = sensorGet(sensor, 'photon rate');
                
        % Compute background adaptation parameters
        bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
        
        initialState = riekeAdaptSteadyState(bgR, p, sz);
        adaptedData  = riekeAdaptTemporal(pRate, initialState);
                       
    otherwise
        error('unknown adaptation type');
end



end
