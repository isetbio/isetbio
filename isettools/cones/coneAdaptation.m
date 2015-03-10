function absorptions = coneAdaptation(typeAdapt,volts,sensor,species)
%Cone adaptation
%
%   absorptions = coneAdaptation(typeAdapt,volts,sensor,species)
%
% Implement adaptation models to produce the cone volts.  The properties of
% the sensor and the type of species are additional parameters.  Not sure
% they should be separated out this way.
%
% The adaptation serves two main goals.  First, the mean level is
% determined as a reference so that increments and decrements are
% meaningful. This is needed so that on- and off-cells can work.  When we
% have a cone absorption that is above the mean, this should cause an
% off-cell to reduce its firing rate and on-cell to increase.  Conversely,
% cone absorptions below the mean cause an on-cell to decrease and an
% off-cell to increase.  By setting a 0 level at the mean, the cones signal
% to the RGCs.  (This text sucks.  Make it better).
%
% See the discussion in Rodieck p. 177 and thereabouts for ideas about cone
% signals.  The physiological question is what is the meaning of the offset
% or zero mean?  Rodieck describes the effect of light as setting the mean
% transmitter release.  In the dark, there is a relatively large dynamic
% range. As the light is on steadily, the release decreases and the range
% available for another flash also decreases.   If you darken from the
% mean, however, the rate can increase.
%
% The way in which the mean is set must depend on a combination of the
% photoisomerizations and the recycling rate, probably through some
% equilibrium equation.  The recycling rate depends on the isomerization
% rate and they set an equilibrium that is slower and slower as the mean
% background gets brighter.  This is how the value of the offset gets set.
%
% In addition to Rodieck's discussion, there are famous papers by Boynton
% and others at SRI expressing such a model based on cone ERPs, and
% probably Rushton and others (van Norren?).
%
% The other issue is the total gain.  The cones can only modulate over a
% dynamic range around the mean.  So, we set the gain as well to keep the
% modulation within some range, such as +/- 50 mV of the mean.
%
% The third issue is whether all the cones are the same, or there is
% space-variant adaptation.  That is for the future.
%
% Inputs
%  typeAdapt:  The adaptation model
%   typeAdapt= 0 - no adaptation
%   typeAdapt= 1 - a single gain map for the whole cone array
%   typeAdapt= 2 - one gain map computed for each cone class
%   typeAdapt= 3 - cone by cone adaptation (NYI)
%
%  volts:      These are a 3D array of coneMosaic x nTimeSamples.
%  sensor:     This is returned with the mean voltage.
%  species:    Normally 'human', but we are experimenting with 'mouse'.
%
% Returns
%  absorptions: A structure containing the adapted data used by the RGC
%  calculations.  Information is added to the structure in the calling
%  routine (coneAbsorptionsTS).
%   .data = the adapted volts
%   .sensor = save the sensor here
%   .unadapteddata = volts;
%   .gain = the effective adaptation gain
%
%
% ALERT:  The amount of adaptation is hard coded at this moment.  It should
% be a parameter somewhere.  Worse yet, the amount of adaptation differs
% between types 1 and 2. Fix this.
%
% In the future, we will look at the time series and have a time-varying
% adaptation function.
%
% See also: coneAbsorptionsTS, rgc_Harmonic
%
% Examples:
%
% (c) Stanford VISTA lab, 2010

%% Why are there 4 adaptation terms?
%  When one is black, it shouldn't become adapated
warning('Use coneAdapt.  This routine (coneAdaptation) is now out of date.')
warning('RGC code still uses this and should be updated')

if notDefined('typeAdapt'), typeAdapt = 2; end
if notDefined('volts'), error('volts are required'); end
if notDefined('sensor'), error('sensor is required'); end
if notDefined('species'), species = 'human'; end

if typeAdapt == 0, adaptedData = volts; gain = 1;
else
    
    % Comments needed for adaptation
    switch lower(species)
        case 'human'
            [adaptedData,gain] = rgcHumanAdapt(volts,sensor,typeAdapt);
        case 'mouse'
            [adaptedData,gain] = rgcMouseAdapt(volts,sensor,typeAdapt);
        otherwise, error('Unknown species %s\n',species);
    end
end

% Set the zero level as the median.  See discussion in the header.
offset      = median(adaptedData(:));
adaptedData = adaptedData - offset;

fprintf('Mean adapted cone response (volts): %f\n',mean(adaptedData(:)))

% The absorptions structure holds adapated data, volts, and other
% structures. This was EC's idea.  Could be simplified and improved.
absorptions.data      = adaptedData;  % Millivolts, both positive and negative
absorptions.typeAdapt = typeAdapt;
absorptions.adaptGain = gain;         % Sometimes a number, sometimes a vector
absorptions.offset    = offset;       % How we set the mean to 0
absorptions.unadapted = volts;        % Millivolts, but only positive.
absorptions.sensor    = sensor;

return

%---------
function [adaptedData, gainMap] = rgcHumanAdapt(volts,sensor,typeAdapt)
%Cone adaptation functions for human
%
% Not sure why it differs from mouse ... but there is some issue about
% wavelength in the mouse.
%
% Two adpatation calculations are implemented.  One gain for all cones, or
% separate amounts for each one
%
% (c) Stanford VISTA lab, 2010

vSwing = pixelGet(sensorGet(sensor,'pixel'),'voltageSwing');

switch typeAdapt
    case 1 % Same gain for all cone types and positions
        
        % Extract the sensor voltages (all of them)
        
        % Set the gain so that the returned mean voltage is 1/2 of voltage
        % swing.  This adaptation has not dynamic, sigh.
        gainMap = (0.8*vSwing) / mean(volts(:));
        adaptedData = gainMap*volts;
        
    case 2 % Method 2 : different gain for each cone type
        % Note : this will help visually, but does it help the SVM?
        % Note : even with no noise sources (even shot noise), the output
        % is not quite "flat", there are still differences, probably because
        % the center is clearer that the edges so the mean doesn't quite
        % work. See figure about this :
        % filename =fullfile(synapseRootPath,'vrn','rgc','data','human_adapt_example.fig');
        % open(filename)
        
        % get cone mosaic
        %     coneMap = sensorGet(sensor,'pattern');
        %     coneMap = coneMap(:); % flatten out to line format
        %
        % this should be 4 : black, S, M, L
        nConeTypes = size(sensorGet(sensor,'filterSpectra'),2);
        
        % The volts come in as multiple samples.  We compute the mean
        sensor = sensorSet(sensor,'volts',mean(volts,3));
        filterLetters = sensorGet(sensor,'filterColorLetters');
        
        % Set a desired level for the return?
        for ii=1:nConeTypes
            v = sensorGet(sensor,'volts',ii);
            m = mean(v(:));
            if m>0 && ~strcmpi(filterLetters(ii),'k')
                % Maybe this should be half of the voltage swing?
                gainMap(ii) = (vSwing/2)/mean(v(:));
            else
                % Black pixel, don't scale.  We should get rid of the black
                % somewhere later in the code
                gainMap(ii) = 1;
            end
        end
        
        % Figure out which cone type is in which position
        [CFAnumbers,mp] = sensorImageColorArray(sensorDetermineCFA(sensor));
        % figure; image(CFAnumbers); colormap(mp)
        
        % determine which value cone type as specified in filterLetters
        % corresponds to which index in CFAnumbers
        inds = sensorImageColorArray(filterLetters);
        
        % Scale each of the cone types by its own gain, keeping the data
        % in place
        adaptedData = volts;
        for ii=1:nConeTypes
            g = ones(size(volts(:,:,1)));
            g(CFAnumbers == inds(ii)) = gainMap(ii);
            g = repmat(g,[1,1,size(adaptedData,3)]);
            adaptedData = adaptedData .* g;
        end
        
    case 3 % Method 3 : adapt each cone independently
        error('Adapt type 3 (independent adaption for each cone) NYI')
        
end

return

% -----
function [adaptedData, gainMap] = rgcMouseAdapt(volts,sensor,typeAdapt)
% Implement cone adaptation for the mouse retina
%
% Code needs more thought and documentation - first draft by EC for class
% project.
%
% The different cones get different quantities of photons, because the lens
% transmittance is different over different wavelength (cuts out some UV).
% In humans, we don't have this problem : transmittance is 1 at all
% wavelengths.
%
% (c) Stanford VISTA lab, 2010

if typeAdapt == 1
    % Absorption gain = 1/number of photons absorbed for equal wavelength illumination
    %                 = 1/integral(filterSpectrum * transmittance)
    % We divide by the mean luminance mLum.
    % Get parameters : wavelengths, spectra of cones/filters, transmittance
    % of lens
    wave = sensorGet(sensor,'wave');
    filterSpectra = sensorGet(sensor,'filterSpectra');
    transmittance = opticsGet(oiGet(oi,'optics'),'transmittance');
    
    % compute gain for each type of cone
    filters = filterSpectra .* repmat(transmittance,1,size(filterSpectra,2));
    integ = sum(filters,1);
    adaptGain = 1./integ;
    
    % apply to the cone map of the retina
    coneMap = sensorGet(sensor,'pattern');
    gainMap = zeros(size(coneMap));
    for coneType = 2:length(adaptGain)
        gainMap = gainMap + (coneMap == coneType)*adaptGain(coneType);
    end
    % apply gain to absorptions
    volts = repmat(gainMap,[1,1,size(volts,3)]) .* volts;
    % => gain too small to compensate!
    %
elseif typeAdapt == 2
    % Method 2--------------
    % Compute mean of each conetype, and add gain to equate means of the two types.
    %
    % Get parameters of cone pattern
    coneMap = sensorGet(sensor,'pattern');
    sz = size(coneMap);
    column = coneMap(:,1);
    numRowsM = length(column(column==column(1)));
    numRowsUV = length(column(column==column(end)));
    numRowsMixed = sz(1) - numRowsM - numRowsUV;
    % get means for each cone type
    meanM = mean(mean(mean(volts(1:numRowsM, :, :))));
    firstRowUV = numRowsM+numRowsMixed+1;
    meanUV = mean(mean(mean(volts(firstRowUV:sz(1), :, :))));
    % Gain : equate all means to M-cone mean
    gainMap = zeros(sz);
    gainMap(1:numRowsM, :) = 1; % no gain for M cones
    gainMap(firstRowUV:sz(1),:) = meanM/meanUV; % gain for UV cones
    % gain vector for mixed cones
    % using means of cone types
    %gainMixed = meanM ./mean(mean(volts((numRowsM+1):(firstRowUV-1),:,:),2),3);
    % using interp of gains
    % gainMixed = interp1([0,numRowsMixed+1]', [1, meanM/meanUV]', [1:numRowsMixed]', 'linear');
    % using interp of means, then inverse
    gainMixed = meanM ./ interp1([0,numRowsMixed+1]', [meanM, meanUV]', (1:numRowsMixed)', 'linear');
    gainMap((numRowsM+1):(firstRowUV-1), :) = repmat(gainMixed, 1, sz(2));
    % Apply gain map to absorption data
    adaptedData = volts .* repmat(gainMap,[1,1,size(volts,3)]);
else
    error('Only two types of adaptation for the mouse');
end

return
