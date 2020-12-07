function [lmsFilters, meanCurrent] = linearFilters(os, cMosaic, varargin)
% Returns the photocurrent impulse response for a single absorption
%
% Syntax:
%	[lmsFilters, meanCurrent] = linearFilters(os, cMosaic)
%
% Description:
%    The LMS impulse response functions calculated here model the cone
%    photocurrent response to brief or low contrasts with respect to a
%    steady background. These experimental conditions (steady backgrounds, 
%    modest contrasts) are often found in psychophysical or physiological
%    experiments.
%
%    The impulse response function is derived from Rieke's biophysical
%    model. It depends on the mean absorption rate, and exhibits adaptation
%    behavior. We take the difference in the response to a constant
%    stimulus and one with a small delta function increment.
%
%    See osLinear.osCompute for how the impulse response function and mean
%    currents are used.
%
%    There are different parameters for foveal and peripheral functions, as
%    implemented by osBioPhys.
%
%    The LMS filters (impulse response functions) are stored here at a
%    particular time step (os.timeStep), which is typically 1 ms, but could
%    be shorter. When it is used in osLinear.osCompute, the filters are
%    resampled to the time base of the cone absorptions.
% 
% Input:
%    os                         - A linear outer segment object
%    cMosaic                    - The parent object of the outersegment
%
% Output:
%    lmsFilters                 - The impulse responses to a single photon
%                                 added to the background; these are stored
%                                 in lmsConeFilter in the os object.
%    meanCurrent                - The steady state caused by the mean
%                                 absorption rate. This value is used in
%                                 osCompute().
%
% Optional key/value pairs:
%    'absorptionsInXWFormat'    - Matrix. Pass input absorptions directly
%                                 in XW format, rather than getting them
%                                 out of the cMosaic object (default
%                                 empty). This is passed into routine
%                                 coneMeanIsomerizations.
%    'eccentricity'             - Value. Eccentricity in degrees to pass on
%                                 to the osBioPhys object when computing
%                                 the linear impulse response. See
%                                 osBioPhys for description of its meaning.
%                                 Default 15.
%
% See Also:
%    coneMeanIsomerizations, osBioPhys, v_osBioPhys, t_coneMosaicFoveal, 
%    t_osLinearize, s_osLinearFilters

%    11/xx/16  JRG/BW  (c) isetbio team
%    08/06/17  dhb     Cleanup comments. 
%                      Explicitly set cone type in the pattern on each
%                      loop, so that that things will work right if we
%                      every put in explicit dynamics by cone type.
%	 02/12/18  jnm     Formatting
%    10/13/2020 npc    Interpolating IR and mean current based on
%                      eccentricity, from data on 0 degs and 10 degs

    %% parse input parameters
    p = inputParser; p.KeepUnmatched = true;
    p.addRequired('os', @(x) isa(x, 'outerSegment'));
    p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic')); 
    p.addParameter('absorptionsInXWFormat', [], @isnumeric);
    p.addParameter('eccentricity', 15, @isnumeric);
    p.parse(os, cMosaic, varargin{:})
    eccentricityDegs = p.Results.eccentricity; 

    %% Get mean isomerization rate, in R*/sec
    meanRate = coneMeanIsomerizations(cMosaic, 'absorptionsInXWFormat', ...
        p.Results.absorptionsInXWFormat);

    % Compute filters and mean current for 0 degs eccentricity
    [lmsFiltersF, meanCurrentF] = biophysicalLinearFilters(0.00, meanRate, os);
    
    % Compute filters and mean current for 15 degs eccentricity
    [lmsFiltersP, meanCurrentP] = biophysicalLinearFilters(15, meanRate, os);
    
    % Interpolate filters and mean currents based on ecc from 0 degs up to osBioPhys.fovealPeripheralCutoffDegs
    % For ecc > osBioPhys.fovealPeripheralCutoffDegs degs, we use the peripheral data
    linearMixingFactor = 1-min([1 eccentricityDegs/osBioPhys.fovealPeripheralCutoffDegs]);
    
    % Compute density-based intepolation factors across eccentricities
    obj = WatsonRGCModel();
    eccDegs = logspace(log10(0.1), log10(osBioPhys.fovealPeripheralCutoffDegs), 64);
    eccUnits = 'deg';
    densityUnits = 'deg^2';
    [~, coneRFDensity] = obj.coneRFSpacingAndDensityAlongMeridian(eccDegs, 'temporal meridian', eccUnits, densityUnits);
    densityBasedMixingFactor  = (coneRFDensity-min(coneRFDensity)) / (max(coneRFDensity)-min(coneRFDensity));
    
    % Density-based factor at the mosaic's ecc
    [~,idx] = min(abs(eccDegs-eccentricityDegs));
    densityBasedMixingFactor = densityBasedMixingFactor(idx);
    
    % Interpolate based on density (non-linear) mixing factor
    mixingFactor = densityBasedMixingFactor;
    lmsFilters  = lmsFiltersF  * (mixingFactor) + (1-mixingFactor) * lmsFiltersP;
    meanCurrent = meanCurrentF * (mixingFactor) + (1-mixingFactor) * meanCurrentP;
end

function [lmsFilters, meanCurrent] = biophysicalLinearFilters(eccentricityDegs, meanRate, os)

    %% Setup parameters
    %
    % We will get the linear impuluse response by explicitly passing in a delta
    % function to the osBioPhys object and extracting what happens. We need to
    % set some actual parameters to make this happen.
    %
    % Parameters
    timeStep = os.timeStep;                % time step (should be < 1 ms)
    nSamples = round(0.8 / timeStep) + 1;  % 0.8 total sec
    flashIntensity = 1;                    % 1 photon above the background mean
    warmupTime = round(0.4 / timeStep);    % Warm up period is 0.4 sec

    % Where we store the filters
    os.lmsConeFilter = zeros(nSamples-warmupTime + 1, length(meanRate));
    meanCurrent = zeros(1, 3);

    %% Generate cone mosaic with an outerSegment based on the biophysical model 
    %
    % We turn off the noise and use the biophysical coneMosaic model to
    % calculate an impulse response. We set up a mosaic with a single L cone as
    % a placeholder, but this gets set to the different types in the loop below
    osCM = osBioPhys('eccentricity', eccentricityDegs);
    % Run it without noise
    osCM.set('noise flag', 'none');
    % Single cone mosaic, L cone as placeholder
    cm = coneMosaic('os', osCM, 'pattern', 2);
    cm.integrationTime = timeStep;
    cm.os.timeStep = timeStep;

    %% For each of the cone types ...
    assert(length(meanRate) == 3, ...
        'We only know how to deal with three cone types');
    for meanInd = 1:length(meanRate)

        % Get the isomerization rate (R*) in each time step R*
        meanIntens  = meanRate(meanInd) * timeStep;  

        % Create a constant stimulus at this rate
        stimulus = meanIntens*ones(nSamples, 1);

        % Compute outer segment current for the constant stimulus The loop goes
        % from 1 to 3 for L, M and S and these cone types are indexed in the
        % pattern by 2, 3, 4. This decision was baked in very early in the
        % ISET/ISETBio design and is hard to back out of now.
        cm.pattern = meanInd + 1;
        cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
        cm.computeCurrent('bgR', meanRate(meanInd));
        currentConstant = squeeze(cm.current);
        % vcNewGraphWin;
        % plot(timeStep * (1:nSamples), currentConstant);

        % Add a single photon (impulse) to the background one time step after
        % the warmup period
        stimulus(warmupTime+1) = meanIntens + flashIntensity;

        % Compute outer segment currents with biophysical model with the
        % impulse stimulus
        cm.absorptions  = reshape(stimulus, [1 1 nSamples]);
        cm.computeCurrent('bgR', meanRate(meanInd));
        currentImpulse = squeeze(cm.current);
        % hold on; plot(timeStep*(1:nSamples), currentImpulse);

        % Store the impulse response. We divide by flashIntensity in for
        % completeness, but it is 1 so really, no need.
        %
        % Impulse response is stored so that it starts at 0 and goes positive
        os.lmsConeFilter(:, meanInd) = ...
            ((currentImpulse((warmupTime:end) - 1)) ...
            - currentConstant((warmupTime:end) - 1)) / flashIntensity;
        % vcNewGraphWin; 
        % plot(stimulus - meanIntens);
        % hold on; 
        % plot(currentImpulse-currentConstant);
        % grid on

        % We are a tiny bit worried about the edges; so we set a few steps
        % before the very end to the mean current from the constant stimulus
        % background. Edges are always a bitch.
        meanCurrent(meanInd) = currentConstant(end-10);
    end

    %% Assign the filters as output and return
    lmsFilters = os.lmsConeFilter;

end