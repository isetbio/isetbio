function varargout = v_wvfDiffractionPSF(varargin)
%
% Checks that diffraction-limited PSFs are correct.
%
% Compares  monochromatic PSFs computed from Zernike coefficients in this
% toolbox with those in PTB and ISET. The curves/points in each lineplot
% figure should overlay.
%
% At the end, we calculate a slice through the PSF for each wavelength.
% Illustrates how to create lines with appropriate spectral colors.
%
% (c) Wavefront Toolbox Team, 2012
    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Initialize
close all; ieInit;

%% Some informative text
UnitTest.validationRecord('SIMPLE_MESSAGE', 'Check diffraction limited PSFs.');

%% Compare pointspread function in wvf with psf in Psych Toolbox

% When the Zernike coefficients are all zero, the wvfComputePSF code should
% return the diffraction limited PSF.  We test whether this works by
% comparing to the diffraction limited PSF implemented in the PTB routine
% AiryPattern.

% Set up default wvf parameters for the calculation 
wvf0 = wvfCreate;

% Specify the pupil size for the calculation
calcPupilMM = 3;
wvf0 = wvfSet(wvf0,'calc pupil size',calcPupilMM);

% Plotting ranges for MM, UM, and Minutes of angle
maxMM = 1;
maxUM = 20;
maxMIN = 2;

% Which wavelength to plot
wList = wvfGet(wvf0,'calc wave');

%% Calculate the PSF, normalized to peak of 1

% This function computes the PSF by first computing the pupil function.  In
% the default wvf object, the Zernicke coefficients match diffraction.
wvf0 = wvfComputePSF(wvf0);

% Make sure psf computed this way (with zcoeffs zeroed) matches
% what is returned by our internal get of diffraction limited psf.
psf = wvfGet(wvf0,'psf');

% We make it easy to simply calculate the diffraction-limited psf of the
% current structure this way.  Here we make sure that there is no
% difference.
diffpsf = wvfGet(wvf0,'diffraction psf');
UnitTest.assertIsZero(max(abs(psf(:)-diffpsf(:))),'Internal computation of diffraction limited psf',0);

% Verify that the calculated and measured wavelengths are the same
calcWavelength = wvfGet(wvf0,'wavelength');
measWavelength = wvfGet(wvf0,'measured wavelength');
UnitTest.assertIsZero(max(abs(measWavelength(:)-calcWavelength(:))),'Measured and calculation wavelengths compare',0);

%% Plots 

% Make a graph of the PSF within maxUM of center
wvfPlot(wvf0,'2dpsf space','um',wList,maxUM);

% Make a graph of the PSF within 2 arc min
wvfPlot(wvf0,'2dpsf angle','min',wList,maxMIN);

%% Plot the middle row of the psf, scaled to peak of 1
wvfPlot(wvf0,'1d psf angle normalized','min',wList,maxMIN);
hold on

% Get parameters needed for plotting comparisons with PTB, below
arcminutes       = wvfGet(wvf0,'psf angular samples','min',wList);
arcminpersample  = wvfGet(wvf0,'ref psf sample interval');
arcminpersample1 = wvfGet(wvf0,'psf arcmin per sample',wList);
arcminpersample2 = wvfGet(wvf0,'psf angle per sample',[],wList);
if (arcminpersample1 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end
if (arcminpersample2 ~= arcminpersample1)
    error('Default units of get on ''psfanglepersample'' unexpectedly changed');
end
index = find(abs(arcminutes) < 2);
radians = (pi/180)*(arcminutes/60);

% Compare to what we get from PTB AiryPattern function -- should match
ptbPSF = AiryPattern(radians,calcPupilMM ,calcWavelength);
plot(arcminutes(index),ptbPSF(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',calcPupilMM,calcWavelength));
UnitTest.validationData('ptbPSF', ptbPSF);

%% Do the same thing using isetbio functions
thisWave = 550;
oi = oiCreate('diffraction limited');
optics = oiGet(oi,'optics');
fLength = 0.017;              % Human focal length is about 17 mm
fNumber = 17/calcPupilMM;     % Set f-number which fixes pupil diameter

optics = opticsSet(optics,'flength',fLength);   % Roughly human
optics = opticsSet(optics,'fnumber',fNumber);   % Roughly human
oi = oiSet(oi,'optics',optics);
uData = oiPlot(oi,'psf',[],thisWave);
set(gca,'xlim',[-10 10],'ylim',[-10 10]);
UnitTest.validationData('oi', oi);

%% Now, compare all three
[r,c] = size(uData.x);
mid = ceil(r/2);
psfMid = uData.psf(mid,:);
posMM = uData.x(mid,:)/1000;               % Microns to mm
posMinutes = 60*(180/pi)*(atan2(posMM,opticsGet(optics,'flength','mm')));

g = wvfPlot(wvf0,'1d psf angle normalized','min',wList,maxMIN);
hold on
plot(posMinutes,psfMid/max(psfMid(:)),'ko')
hold on
plot(arcminutes(index),ptbPSF(index),'b','LineWidth',2);
xlabel('Arc min')
set(gca,'xlim',[-2 2])
grid on
legend('WVF','ISETBIO','PTB');

UnitTest.validationData('wvf0', wvf0);

%% Repeat the PSF calculation with a wavelength offset

% This section checks that if we add an explicit observer focus correction,
% in this case the amount needed to correct for chromatic aberration, we
% get the same result.  It is a pretty small test of the function
% wvfLCAFromWavelengthDifference relative to the measured wavelength

% Copy the wavefront structure
wvf1 = wvf0;
wvf17 = wvf0;

% Let's work at this very short wavelength
wList = 400;
wvf1 = wvfSet(wvf1,'calc wave',wList);

% This is the chromatic aberration relative to the measured wavelength
lcaDiopters = wvfLCAFromWavelengthDifference(wvfGet(wvf1,'measured wl'),wList);

%  We set the parameter as if the measurement has this correction.
wvf1 = wvfSet(wvf1,'calc observer focus correction',lcaDiopters);
wvf1 = wvfComputePSF(wvf1);
w = wvfGet(wvf1,'calc wave');
pupilSize = wvfGet(wvf1,'calc pupil size','mm');

% Get diffractoin limited PSF with both measured and calc wavelength set
% the same and to the short wavelength. This should give us the same
% diffraction limited answer at the short wavelength as the other way of
% computing above (and is the more naturally expressive way to get this
% through the wavefront method, if that is all you are trying to do.)
%
% Fuss here with the wvf structure to up spatial sampling density on this
% one so we get a smooth comparison with the others.  The spatial sampling
% in the spatial domain depends on the measurement wavelength, so we don't
% have the psf at exactly the same spatial points here as in the others.
% By upping the density we can look in the plot and see the smooth curve
% pass through the more coarsely computed points.
%
% In the plot created below, the result of this one is the thin green line, which 
% you can examine and see passing through the sample points for the red and blue
% curves, which are more coarser.
wvf17Samples = 1001;
wvf17 = wvfSet(wvf17,'measured wl',wList);
wvf17 = wvfSet(wvf17,'calc wave',wList);
wvf17 = wvfSet(wvf17,'number spatial samples',wvf17Samples);
wvf17 = wvfSet(wvf17,'ref psf sample interval',wvfGet(wvf17,'ref psf sample interval')/4);
wvf17 = wvfComputePSF(wvf17);

% There should be no difference (again) because we corrected for the
% chromatic aberration.  As noted above, the thin green curve is smoother
% and deviates in between the coarser samples for the red and blue curves,
% but you can see that the function being computed is in good agreement.
wvfPlot(wvf1,'1d psf angle normalized','min',w,maxMIN);
hold on
[~,h] = wvfPlot(wvf17,'1d psf angle normalized','min',w,maxMIN,'no window');
set(h,'Color','g','LineWidth',1);
ptbPSF = AiryPattern(radians,pupilSize,w);
plot(arcminutes(index),ptbPSF(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalize PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,w));

% Save unit test data
UnitTest.validationData('wvf1', wvf1);
UnitTest.validationData('wvf17', wvf1);
UnitTest.validationData('ptbPSF1', ptbPSF);

% PSF angular sampling should be the same across wavelengths
arcminpersample2 = wvfGet(wvf1,'psf angle per sample','min',w);
if (arcminpersample2 ~= arcminpersample)
    error('PSF sampling not constant across wavelengths');
end

%% Use a different pupil size at original wavelength

% Copy the original
wvf2  = wvf0;

% Calculate for a larger pupil (less diffraction, therefore.
pupilMM = 7; 
wvf2  = wvfSet(wvf2,'calc pupil diameter',pupilMM);

% Confirm parameters
wvf2  = wvfComputePSF(wvf2);
wList = wvfGet(wvf2,'calc wave');
pupilSize = wvfGet(wvf2,'calc pupil size','mm');

% Compare the PTB and WVF curves
wvfPlot(wvf2,'1d psf angle normalized','min',wList,maxMIN);
ptbPSF = AiryPattern(radians,pupilSize,wList);

hold on
plot(arcminutes(index),ptbPSF(index),'b','LineWidth',2);
xlabel('Arc Minutes');
ylabel('Normalized PSF');
title(sprintf('Diffraction limited, %0.1f mm pupil, %0.f nm',pupilSize,wList));

UnitTest.validationData('wvf2', wvf2);
UnitTest.validationData('ptbPSF2', ptbPSF);

%% Show the PSF slices across wavelengths along with the 'white'

% New copy
wvf3 = wvf0;

% This makes a colormap that is like the spectral colors
pupilMM  = 3.0;
thisWave = 420:10:650;
cmap = squeeze(xyz2srgb(XW2RGBFormat(ieReadSpectra('XYZ',thisWave),length(thisWave),1)));

% We compare many the wavelengths and the average across them (white)
wvf3 = wvfSet(wvf3,'calc wave',thisWave);
wvf3 = wvfSet(wvf3,'calc pupil diameter',pupilMM);
wvf3 = wvfComputePSF(wvf3);

%vcNewGraphWin([],'tall'); 
vcNewGraphWin;
for ii=1:length(thisWave)
    if ii == 1
        [u,pData] = wvfPlot(wvf3,'1d psf space','um',thisWave(1),5*maxMIN,'no window');
        x = u.x; y = u.y/sum(u.y(:));
        set(pData,'color',cmap(ii,:),'LineWidth',1);
    end
    hold on
    [uData, pData] = wvfPlot(wvf3,'1d psf space','um',thisWave(ii),'no window');
    thisY = interp1(uData.x,uData.y,x);
    y = y + thisY;
    set(pData,'color',cmap(ii,:),'LineWidth',1);
end
str = num2str(thisWave');

% Calculate the average
y = y/length(thisWave);
p = plot(x,y,'k:'); set(p,'LineWidth',2);
str(end+1,:) = 'wht';

% Uncomment this to see diffraction, but it is the same as 550nm
% hold on
% ptbPSF = AiryPattern(radians,pupilSize,w);
% plot(arcminutes(index),ptbPSF(index),'k:','LineWidth',1);
% str(end+1,:) = 'dfl';

% Labels
xlabel('Position (um)');
ylabel('Slice through PSF');
set(gca,'xlim',[-20 20])
title(sprintf('DL %0.1f mm pupil (water)',wvfGet(wvf3,'calc pupil diameter')));

end


