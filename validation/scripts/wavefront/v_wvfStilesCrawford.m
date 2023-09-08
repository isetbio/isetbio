function varargout = v_wvfStilesCrawford(varargin)
%
% Valide Stiles-Crawford effect object.
%
% Note from DHB.  I'm not sure here what the effect should look like, so
% at present this simply validates that the code runs witout error and
% that the SCE has an effect on the PSF.  The PSF does get narrower with
% the SCE, which is expected.
%
% 8/19/12  dhb  Updated.
% 8/12/15  dhb  UnitTestToolbox'd.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Initialize ISETBIO
close all;
toleranceFraction = 0.001;

%% Some informative text
UnitTest.validationRecord('SIMPLE_MESSAGE', 'Validate wavefront Stiles-Crawford Effect code.');

%% For plotting limits
maxMIN = 2;
maxMM  = 1;
waveIdx = 1;
theWavelength = 550;

%% Set up wvf structure and sce structure
wvf = wvfCreate;
wvf = wvfSet(wvf,'zcoeffs',[0.2 0.75],{'defocus', 'oblique_astigmatism'});
sceP = sceCreate(theWavelength,'berendschot_data','centered');

%% No Stiles Crawford effect
wvf = wvfSet(wvf,'sce params',[]);
wvf = wvfComputePupilFunction(wvf);
wvf = wvfComputePSF(wvf);
sce1DFig2 = vcNewGraphWin; hold on
wvfPlot(wvf,'1d psf angle','min',[],maxMIN,'no window');

% Get variables to validate
zCoeffs = wvfGet(wvf,'zcoeffs');
theTolerance = mean(zCoeffs(:))*toleranceFraction;
UnitTest.validationData('zCoeffsNoSCE', zCoeffs, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'zCoeffsNoSCE', theTolerance);

pupilFunction = abs(wvfGet(wvf,'pupil function'));
theTolerance = mean(pupilFunction(:))*toleranceFraction;
UnitTest.validationData('pupilFunctionNoSCE', pupilFunction, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'pupilFunctionNoSCE', theTolerance);

PSF = wvfGet(wvf,'PSF');
theTolerance = mean(PSF(:))*toleranceFraction;
UnitTest.validationData('PSFNoSCE', pupilFunction, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'PSFNoSCE', theTolerance);

%% Include the SCE in place
wvf = wvfSet(wvf,'sce params',sceP);
wvf = wvfComputePSF(wvf);
[f,p] = wvfPlot(wvf,'1d psf angle','min',[],maxMIN,'no window');
set(p,'color','b')
hold on

% Get variables to validate
zCoeffs = wvfGet(wvf,'zcoeffs');
theTolerance = mean(zCoeffs(:))*toleranceFraction;
UnitTest.validationData('zCoeffsWithSCE', zCoeffs, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'zCoeffsWithSCE', theTolerance);

pupilFunction = abs(wvfGet(wvf,'pupil function'));
theTolerance = mean(pupilFunction(:))*toleranceFraction;
UnitTest.validationData('pupilFunctionWithSCE', pupilFunction, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'pupilFunctionWithSCE', theTolerance);

PSF = wvfGet(wvf,'PSF');
theTolerance = mean(PSF(:))*toleranceFraction;
UnitTest.validationData('PSFWithSCE', pupilFunction, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'PSFWithSCE', theTolerance);

end

