%% v_wvfPupilFunction
%
% Explore the effect of changing the pupil size in the calculation.
%
% We load the Thibos wavefront data collected at some pupil diameter.
% There are a few options.  We set these data into the zernicke
% coefficients.
%
% We then set the calculated pupil size and compute the expected
% pointspread function.  That function changes as the pupil diameter gets
% smaller.
%
% (BW) (c) Wavefront Toolbox Team, 2014

%% 
s_initISET

%% Load the Thibos data for one of the pupil diameter sizes 
pupilMM = 7.5;   % Could be 6, 4.5, or 3
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);

% Create the wvf parameter structure with the appropriate values
wave = (400:10:700)';
wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('%d-pupil',pupilMM));
wvfP = wvfSet(wvfP,'measured pupil',pupilMM);

%% Calculate the effect of varying the pupil diameter
cPupil = [2,3,4,5,6,7];
for ii=1:sum(cPupil<=pupilMM)
    wvfP = wvfSet(wvfP,'calculated pupil',cPupil(ii));
    wvfP = wvfComputePSF(wvfP);
    wvfPlot(wvfP,'2d psf space','um',550,20)
    title(sprintf('Calculated pupil diameter %.1f mm',cPupil(ii)));
end

%% End