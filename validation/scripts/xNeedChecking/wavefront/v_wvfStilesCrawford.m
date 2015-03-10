%% v_wvfStilesCrawford
%
% Note from DHB.  I'm not sure here what the effect should look like, so
% at present this simply validates that the code runs witout error and
% that the SCE has an effect on the PSF.
%
% 8/19/12  dhb  Updated.

%% Initialize
s_initISET

% For plotting limits
maxMIN = 2;
maxMM  = 1;
waveIdx = 1;
theWavelength = 550;

%% Set up wvf structure and sce structure
wvfParams = wvfCreate;
wvfParams = wvfSet(wvfParams,'zcoeffs',[0.2 0.75],{'defocus', 'oblique_astigmatism'});
sceP = sceCreate(theWavelength,'berendschot_data','centered');

%% No Stiles Crawford effect
wvfParams = wvfSet(wvfParams,'sce params',[]);
wvfParams = wvfComputePSF(wvfParams);
sce1DFig2 = vcNewGraphWin; hold on
wvfPlot(wvfParams,'1d psf angle','min',[],maxMIN,'no window');

%%  Include the SCE in place
wvfParams = wvfSet(wvfParams,'sce params',sceP);
wvfParams = wvfComputePSF(wvfParams);
[f,p] = wvfPlot(wvfParams,'1d psf angle','min',[],maxMIN,'no window');
set(p,'color','b')
hold on

%% End

