%% PSF 2 Wernike Coefficients
%
% Search for the Wernike polynomial coefficients that produce a target PSF.
% 
% It is possible there is an analytical solution to this.  I was too lazy
% to figure that out.  I will check with someone over coffee.
%
% Note this good Watson tutorial.  But it has no units, wavelength,
% and such.
%
%  http://jov.arvojournals.org/article.aspx?articleid=2213266
%
% BW, Vistasoft team, 2018

%% This works well in ISETCAM but not in ISETBio
%
% The wvf structure and sets/gets in ISETBio are simplified so that I can
% understand them.  In ISETBio, there are measured, calc, as well as
% different spatial sampling domains. Somewhere 

%% Create a wavefront object with some coefficients

% Create a wvf object
wave = 550;
wvf = wvfCreate('wave',wave);
wvf = wvfSet(wvf,'zcoeffs',.2,'defocus');
wvf = wvfSet(wvf,'zcoeffs',0,'vertical_astigmatism');

% This is needed to make the pupilPlaneSizeMM agree with ISETCam
% But I don't understand it.
wvf = wvfSet(wvf,'sample interval domain','pupil'); 
wvf = wvfSet(wvf,'z pupil diameter',8); 
% wvfGet(wvf,'ref pupil plane size')
% wvf = wvfSet(wvf,'ref pupil plane size',16.2120);
wvf = wvfSet(wvf,'ref pupil plane size',50);
% wvf = wvfSet(wvf,'ref pupil plane size',25);

wvf = wvfComputePSF(wvf);

% This plot is not the same in ISETBio and ISETCam
% wvfPlot(wvf,'image pupil phase','mm');
wvfPlot(wvf,'image psf space','um')
% wvfPlot(wvf,'image pupil amp','mm')

%% Get the parameters we need for the search

thisWaveUM  = wvfGet(wvf,'wave','um');
thisWaveNM  = wvfGet(wvf,'wave','nm');
pupilSizeMM = wvfGet(wvf,'pupil diameter','mm');
zpupilDiameterMM = wvfGet(wvf,'z pupil diameter');

pupilPlaneSizeMM = wvfGet(wvf,'pupil plane size','mm',thisWaveNM);
nPixels = wvfGet(wvf,'spatial samples');
wvf     = wvfComputePSF(wvf);

% These are the spatial samples and psf value.
% When we get a psf from another source, we should interpolate the
% values to these spatial samples.
samp      = wvfGet(wvf, 'psf spatial samples', 'um', wave);
psfTarget = wvfGet(wvf,'psf');   % The spatial sample positions

%{
% This example interpolates the data from t_gullstranEyeTrace into
% the spatial samples required for searching.  We should be able to
% use this logic
%
% I ran the t_gullstranEyeTrace code to get the oi
% Then ...
 illuminance = oiGet(oi,'illuminance');
 s = oiGet(oi,'spatial support','um'); 
 tmp = interp2(s(1,:,1),s(:,1,2),illuminance,samp,samp(:),'cubic',0);
 vcNewGraphWin; mesh(samp,samp,tmp)
 
 % Set this
 psfTarget = tmp/sum(tmp(:));
 %  now run.
%}
f = @(x) psf2zcoeff(x,psfTarget,pupilSizeMM,zpupilDiameterMM,pupilPlaneSizeMM,thisWaveUM, nPixels);

% I should to figure out how to set the tolerances.  Default is 1e-4
zcoeffs = wvfGet(wvf,'zcoeffs');

% I am searching over the first 6 coefficients.  This includes defocus.
% Could do more, I suppose.  Also, the first coefficient ('piston') has no
% impact on the PSF.  So the search always forces that to 0.
nCoeffs = 6;
zcoeffs(1:nCoeffs)
x0 = zeros(size(zcoeffs(1:nCoeffs)));
options = optimset('PlotFcns',@optimplotfval);

x = fminsearch(f,x0,options);

% Piston comes back as an arbitrary value because the error function
% ignores it. We force it to zero here.
x(1) = 0;  

%% Compare the values
fprintf('Estimated zcoeffs\n');
disp(x)

% Compare the PSFs to make sure we have a match
wvf2 = wvfSet(wvf,'zcoeffs',x);
wvf2 = wvfComputePSF(wvf2);
wvfPlot(wvf2,'image psf space','um')
title('Estimated PSF');

vcNewGraphWin;
imagesc(samp,samp,psfTarget)
title('Target PSF'); axis image; colormap(hot)
grid on

%
psf = wvfGet(wvf2,'psf',wave);
vcNewGraphWin;
plot(psfTarget(1:5:end),psf(1:5:end),'.');
axis equal; identityLine;
grid on;

%%
fprintf('True zcoeffs\n');
disp(zcoeffs(1:nCoeffs))

%% Show the pupil phase functions

vcNewGraphWin([],'tall');
subplot(2,1,1), wvfPlot(wvf,'image pupil phase','mm',wave,'no window')

wvf2 = wvfSet(wvf,'zcoeffs',x);
wvf2     = wvfComputePSF(wvf2);
subplot(2,1,2), wvfPlot(wvf2,'image pupil phase','mm',wave,'no window')

%%