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
    
    deltaPhaseDegs = 2;
    phaseDegs = 0:deltaPhaseDegs:(360-deltaPhaseDegs);
    phaseRadians = phaseDegs/180*pi;
    
    deltaAmplitude = 0.01;
    amplitudes = 0.8:deltaAmplitude:1.5;
    
    fTime = 2.0 * pi * temporalFrequencyHz * temporalSupportSeconds;
    theResiduals = inf(numel(phaseDegs), numel(amplitudes));
   
    % Search over phase
    for iPhase = 1:numel(phaseDegs)
        thePhaseRadians = phaseRadians(iPhase);
        theSine = sin(fTime - thePhaseRadians);

        % Search over amplitude
        for iAmp = 1:numel(amplitudes)
            theAmplitude = amplitudes(iAmp);
            theSinewave = theAmplitude * theSine;
            theDifference = theSinewave(:)-theResponse(:);
            r2 = sum(theDifference(:).^2,1);
            theResiduals(iPhase,iAmp) = r2;
        end
    end
    
    [~,idx] = min(theResiduals(:));
    [iPhase, iAmp] = ind2sub(size(theResiduals), idx);
    thePhaseDegs = phaseDegs(iPhase);

    fittedParams = [amplitudes(iAmp)*maxAmplitude thePhaseDegs meanResponse];
    
    % Generate high-resolution fitted function
    fTime = 2.0 * pi * temporalFrequencyHz * temporalSupportHR;
    theFittedResponse = meanResponse + fittedParams(1) * sin(fTime - fittedParams(2)/180*pi);
end
