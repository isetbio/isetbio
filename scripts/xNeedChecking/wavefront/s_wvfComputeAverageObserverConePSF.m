% wvfComputeAverageObserverConePSF
%
%  OBSOLETE - NEEDS TO BE REWRITTEN
%
% Compute the cone PSFs of an average observer.  Done by combining measurements
% from a dataset of Zernike coefficients for a number of observers.
%
% This is done with respect to a specified spectral weighting function,
% and focus optimzed for a specified weighting of the different cone classes.
%
% See also: ComputeConePSFTest, wvfComputeConePSF
%   ComputePSFTest, wvfComputePSF, wvfComputePupilFunction,
%   sceGetParamsParams, wvfGetDefocusFromWavelengthDifference
%
% 3/14/12  MDL  Changed most struct field assignments to use wvfSet(). Updated
%               to set fieldSampleRateMMperPixel instead of sizeOfFieldPixels
% 8/29/11  dhb  Wrote it.

%% Clear
error('This function is obsolete.');

clear; close all;

%% Load cone sensitivities, set weighting spectrum.
S = [400 5 61];
wls = SToWls(S);
load T_cones_ss2;
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);
load spd_D65
weightingSpectrum = SplineSpd(S_D65,spd_D65,S);

% Specify datafile for Zernike coefficients
zernikeFile = 'sampleZernikeCoeffs.txt';
measpupilMM = 6;
theZernikeCoeffs = importdata(zernikeFile);

wvfParams0 = wvfCreate;
wvfParams0 = wvfSet(wvfParams0,'measured pupil',measpupilMM);
wvfParams0 = wvfSet(wvfParams0,'calculated pupil',3);
wvfParams0 = wvfSet(wvfParams0,'wave',wls);
wvfParams0 = wvfSet(wvfParams0,'nominalFocusWl',550);
wvfParams0 = wvfSet(wvfParams0,'defocusDiopters',0);
wvfParams0 = wvfSet(wvfParams0,'fieldSampleSize',16.212/201);
wvfParams0 = wvfSet(wvfParams0,'fieldsizemm',16.212);
wvfParams0.T_cones = T_cones;
wvfParams0 = wvfSet(wvfParams0,'weightspectrum',weightingSpectrum);
whichRow = floor(wvfGet(wvfParams0,'npixels')/2) + 1;

plotLimit = 2;
DOSCE = 1;
if (DOSCE)
    wvfParams0 = wvfSet(wvfParams0,'sceparams',sceCreate(wls,'berendschot_datac'));
else
    wvfParams0 = wvfSet(wvfParams0,'sceparams',sceCreate(wls,'none'));
end
CIRCULARLYAVERAGE = 1;
wvfParams0.coneWeights = [1 1 0];
wvfParams0.criterionFraction = 0.9;

% Read coefficients and optimze PSF for each observer
for i = 1:size(theZernikeCoeffs,2)
    wvfParams = wvfParams0;
    wvfParams.zcoeffs = theZernikeCoeffs(:,i);
    wvfParams = wvfComputeOptimizedConePSF(wvfParams);
    arcminperpixel(i) = wvfParams.arcminperpix;
    defocusDiopters(i) = wvfParams.defocusDiopters;
    lpsfo(:,:,i) = psfCenter(wvfParams.conepsf(:,:,1));
    mpsfo(:,:,i) = psfCenter(wvfParams.conepsf(:,:,2));
    spsfo(:,:,i) = psfCenter(wvfParams.conepsf(:,:,3));
end

% Get optimized diffrac limited PSF
wvfParams = wvfParams0;
wvfParams.zcoeffs = zeros(61,1);
wvfParams = wvfComputeOptimizedConePSF(wvfParams);
arcminperpixeld = wvfParams.arcminperpix;
defocusDioptersd = wvfParams.defocusDiopters;
lpsfd = psfCenter(wvfParams.conepsf(:,:,1));
mpsfd = psfCenter(wvfParams.conepsf(:,:,2));
spsfd = psfCenter(wvfParams.conepsf(:,:,3));

% Get average LMS PSFs
avglpsfo = psfAverageMultiple(lpsfo);
avgmpsfo = psfAverageMultiple(mpsfo);
avgspsfo = psfAverageMultiple(spsfo);
if (CIRCULARLYAVERAGE)
    avglpsfo = psfCircularlyAverage(avglpsfo);
    avgmpsfo = psfCircularlyAverage(avgmpsfo);
    avgspsfo = psfCircularlyAverage(avgspsfo);
    lpsfd = psfCircularlyAverage(lpsfd);
    mpsfd = psfCircularlyAverage(mpsfd);
    spsfd = psfCircularlyAverage(spsfd);
end
onedLPSFo = avglpsfo(whichRow,:);
onedMPSFo = avgmpsfo(whichRow,:);
onedSPSFo = avgspsfo(whichRow,:);
onedLPSFd = lpsfd(whichRow,:);
onedMPSFd = mpsfd(whichRow,:);
onedSPSFd = spsfd(whichRow,:);
maxY = max(max([onedLPSFo(:) onedMPSFo(:) onedSPSFo(:) onedLPSFd(:) onedMPSFd(:) onedSPSFd(:)]));

arcminutes = arcminperpixel(1)*((1:wvfGet(wvfParams0,'npixels'))-whichRow);
index = find(abs(arcminutes) < plotLimit);
figure; clf;
subplot(1,3,1); hold on
plot(arcminutes(index),onedLPSFo(index),'r','LineWidth',4);
plot(arcminutes(index),onedLPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
ylim([0 maxY]);
if (CIRCULARLYAVERAGE)
    title('Circularized L cone PSF');
else
    title('L cone PSF');
end
subplot(1,3,2); hold on
plot(arcminutes(index),onedMPSFo(index),'g','LineWidth',4);
plot(arcminutes(index),onedMPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
ylim([0 maxY]);
if (CIRCULARLYAVERAGE)
    title('Circularized M cone PSF');
else
    title('M cone PSF');
end
subplot(1,3,3); hold on
plot(arcminutes(index),onedSPSFo(index),'b','LineWidth',4);
plot(arcminutes(index),onedSPSFd(index),'k','LineWidth',4);
xlabel('Arc Minutes');
ylabel('PSF');
ylim([0 maxY]);
if (CIRCULARLYAVERAGE)
    title('Circularized S cone PSF');
else
    title('S cone PSF');
end
drawnow;

% Save
if (DOSCE)
    if (CIRCULARLYAVERAGE)
        outfile = sprintf('AverageConePSFCIRC_SCE_%d_%d_%d_%0.1f',wvfParams0.coneWeights(1),wvfParams0.coneWeights(2),wvfParams0.coneWeights(3),wvfParams0.criterionFraction);
    else
        outfile = sprintf('AverageConePSF_SCE_%d_%d_%d_%0.1f',wvfParams0.coneWeights(1),wvfParams0.coneWeights(2),wvfParams0.coneWeights(3),wvfParams0.criterionFraction);
    end
else
    if (CIRCULARLYAVERAGE)
        outfile = sprintf('AverageConePSFCIRC_%d_%d_%d_%0.1f',wvfParams0.coneWeights(1),wvfParams0.coneWeights(2),wvfParams0.coneWeights(3),wvfParams0.criterionFraction);
    else
        outfile = sprintf('AverageConePSF_%d_%d_%d_%0.1f',wvfParams0.coneWeights(1),wvfParams0.coneWeights(2),wvfParams0.coneWeights(3),wvfParams0.criterionFraction);
    end
end
save(outfile,'avglpsfo','avgmpsfo','avgspsfo','lpsfd','mpsfd','spsfd','arcminperpixel','defocusDiopters');





