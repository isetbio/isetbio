%
% RGCMosaicConstructor.helper.fit.sinusoidToResponseTimeSeries
%

function [theFittedResponse, fittedParams] = sinusoidToResponseTimeSeries(...
    temporalSupportSeconds, theResponseTimeSeries, temporalFrequencyHz, temporalSupportHR, varargin)

    p = inputParser;
    p.addParameter('allowOffset', false, @islogical);
    % Execute the parser
    p.parse(varargin{:});
    allowOffset = p.Results.allowOffset;

    if (allowOffset)
        meanResponse = mean(theResponseTimeSeries(:));
    else
        meanResponse = 0;
    end

    theResponseTimeSeries = theResponseTimeSeries - meanResponse;

	maxAmplitude = max(abs(theResponseTimeSeries));
    theResponse = double(theResponseTimeSeries/maxAmplitude);
    
    phaseDegs = 0:2:179;
    phaseRadians = phaseDegs/180*pi;
    
    deltaAmplitude = 0.01;
    amplitudes = 0.8:deltaAmplitude:1.5;
    
    fTime = 2.0 * pi * temporalFrequencyHz * temporalSupportSeconds;
    theResiduals = inf(numel(phaseDegs)*2, numel(amplitudes));
   
    % Search over phase
    for iPhase = 1:numel(phaseDegs)
        thePhaseRadians = phaseRadians(iPhase);
        theSine = sin(fTime - thePhaseRadians);

        % Search over amplitude
        for iAmp = 1:numel(amplitudes)
            theAmplitude = amplitudes(iAmp);
            theSinewave = theAmplitude * theSine;

            theResiduals(iPhase,iAmp) = sum((theSinewave(:)-theResponse(:)).^2);
            theResiduals(iPhase+numel(phaseDegs),iAmp) = sum((-theSinewave(:)-theResponse(:)).^2);
        end
    end
    
    [~,idx] = min(theResiduals(:));
    [iPhase, iAmp] = ind2sub(size(theResiduals), idx);
    if (iPhase > numel(phaseDegs))
        thePhaseDegs = phaseDegs(iPhase-numel(phaseDegs))+180;
    else
        thePhaseDegs = phaseDegs(iPhase);
    end
    fittedParams = [amplitudes(iAmp)*maxAmplitude thePhaseDegs meanResponse];
    
    % Generate high-resolution fitted function
    fTime = 2.0 * pi * temporalFrequencyHz * temporalSupportHR;
    theFittedResponse = meanResponse + fittedParams(1) * sin(fTime - fittedParams(2)/180*pi);
    
    
    debug = ~true;
    if (debug)

        maxR = max([1 max(abs(theFittedResponse(:))) max(abs(theResponse(:) * maxAmplitude))]);

        figure(123); clf;
        plot(temporalSupportSeconds, theResponse * maxAmplitude, 'ks');
        hold on;
        plot(temporalSupportHR, theFittedResponse, 'r-');
        set(gca, 'YLim', maxR*[-1 1]);
        pause
    end 
end
