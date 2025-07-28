%
% RGCMosaicAnalyzer.compute.unwrappedPhaseOfBandPassSTF(theComplexSTF, mode, forceBandPassSolution)
%
function [unwrappedPhase,isBandPass] = unwrappedPhaseOfBandPassSTF(theComplexSTF, mode, forceBandPassSolution)
% Complex STF in, unwrapped phase out. 
%     mode = 'Deg', or 'Rad'
%     forceBandPassSolution: do not detect bandpassiness, force it
%     unwrappedPhase = new unwrapping algorithm bandpass data friendly
%     isBandPass = flag that indicates that the input data was classified as being bandpass.
% 
% Unwrapping bandpass filter phase can be erratic. 
%
% The stock phase unwrap function unwrap(phase) begins its work at the 
% start of the phase array and proceeds to the end. 
%
% With bandpass data, in particular real-world MEASURED bandpass phase data, 
% the beginning of the data set is in the stopband of the filter and noisy. 
% The stock phase unwrapping algorithm therefore bases its initial phase on poor 
% quality data. The noise creates phase offsets that make phase comparisons
% of similar bandpass filter difficult. 
%
% The algorithm works its way out from the middle of the data set to 
% each end rather than starting at the beginning. 
% NOTE: lowpass, highpass, and bandreject data *should* work the same as the 
% stock function but have not been well tested.  Trust but Verify.
% 
% NPC Modified the version of Dick Benson, January 2022
    
    switch mode
        case 'Deg'
            scale = 180/pi;
        case 'Rad'
            scale = 1;
        otherwise
            disp('Specify Deg or Rad')
            return;
    end

    STFmagnitude = abs(squeeze(theComplexSTF));        % magnitude
    STFphaseRadians = angle(squeeze(theComplexSTF));   % phase, not unwrapped.
    unwrappedPhase = unwrap(STFphaseRadians);          % Matlab's unwrap()
    
    N = length(STFphaseRadians);
    if N < 10
       % Not enough data for the new algorithm to be viable. 
       fprintf(2, 'Need at least 10 points to unwrap the STF phase. Using Matlab unwrap\n')
       unwrappedPhase = scale*unwrappedPhase;
       isBandPass = false;
       return
    end
      
    isBandPass = [];
    if (~forceBandPassSolution)
        % Test to see if the data has bandpass characteristics. 
        NP = 3; % 1/2 the data window range that determines M_edges and M_center. 
        
        % Finding the center of the (potentially) bandpass data is the key to 
        % stable results. The Center of Mass calculation works well. 
        centerPointIndex  = round(centerOfMass(STFmagnitude)); 
        
        % If ratio of the magnitude of the in-band complex data (M_center)  
        % to the band edges (M_edges) is greater than a Threshold the data
        % is deemed to be bandpass in nature and the alternative unwrap algorithm is used.
        M_center  = sum(STFmagnitude(  (centerPointIndex-NP):(centerPointIndex+NP) ));
        M_edges   = sum(STFmagnitude(1:NP)) + sum(STFmagnitude( (N-NP+1):N));
        Mag_ratio    = M_center/M_edges;

        % Variations on the above could likely detect highpass data, but this is 
        % left for v2.0.
        Threshold = 20;  % heuristically determined threshold and subject to change. 

        isBandPass = Mag_ratio > Threshold;
    else
        [~,centerPointIndex] = max(STFmagnitude); 
    end

    
    if (forceBandPassSolution || (~isempty(isBandPass)&&isBandPass))
        % Now to unwrap the BPF  phase ... 
        % 1) The middle of the passband is located at the centerPointIndex.
        % 2) Start from centerPointIndex  and work to the right (+freq) with stock unwrap()
        % 3) Flip the data from centerPointIndex to the beginning and unwrap this segment.
        % 4) Reflip the results of step 3 and then append 2&3.
        % 5) Apply scaling for Degrees (x180/pi or Rads (x1)
        phaseRight = unwrap( STFphaseRadians(centerPointIndex:N) );
        phaseLeft  = flip( unwrap( flip( STFphaseRadians(1:centerPointIndex)) ) );
        
        [rows,cols] = size(STFphaseRadians);
        if cols == 1
           % column vector format
           unwrappedPhase = scale*[phaseLeft(1:(centerPointIndex-1)); phaseRight];
        else
           % row vector format
           unwrappedPhase = scale*[phaseLeft(1:(centerPointIndex-1)), phaseRight]; 
        end
        
        fprintf('Unwrapped phase using bandpass solution\n')
    else
        % The data FAILS the test for being bandpass in nature. 
        isBandPass = 0;
        unwrappedPhase = scale * unwrappedPhase; % return stock unwrapped phase.
        fprintf('Unwrapped phase using MATLAB solution\n')
    end
end
    
    
function y = centerOfMass(x)  
  % This function treats the elements of the input vector x as mechanical masses distributed 
  % along a straight line. It computes the location of the balance point.
    
  acc = 0;
  for k = 1:length(x)
      acc = acc + k*x(k);
  end
  y = acc/sum(x);

end