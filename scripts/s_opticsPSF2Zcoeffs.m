%% s_opticsPSF2Zcoeffs
%
% Start with a target PSF and estimate the Zernike coefficients that
% produce that target PSF.
%
% Note this good Watson tutorial.  But it has no units, wavelength,
% and such.
%
%    http://jov.arvojournals.org/article.aspx?articleid=2213266
%
% Current status - 
%  Works for small values, but not for even moderately big ones.
%  The search doesn't get close.  Not sure about which coefficients
%  are better or worse.
%
% BW, Vistasoft team, 2018

%% This works well in ISETCAM but not in ISETBio
%
% The wvf structure and sets/gets in ISETCam are simplified so that I can
% understand them.  In ISETBio, there are measured, calc, as well as
% different spatial sampling domains.

ieInit

%% Create a wavefront object with some coefficients

% Create a wvf object
wave = 550;
wvf = wvfCreate('wave',wave);

% In ISETCam we refer to the z pupil diameter, which means the
% diameter that is attached to the zcoeffs.  In ISETBio the zcoeffs
% refer to the measured pupil diameter.
wvf = wvfSet(wvf,'measured pupil diameter',8);

% Set the defocus coefficient
% [D,V] = meshgrid( 0.3:0.4:0.7, 0.2:.1:.3);
[D,V] = meshgrid( 0.3:0.4:1.0, 0);

%%
pList = [D(:),V(:)];
for ii=1:size(pList,1)
    wvf = wvfSet(wvf,'zcoeffs',pList(ii,1),'defocus');
    % Set the vertical astigmatism
    wvf = wvfSet(wvf,'zcoeffs',pList(ii,2),'vertical_astigmatism');
    
    % This is needed to make the pupilPlaneSizeMM agree with ISETCam
    wvf = wvfSet(wvf,'sample interval domain','pupil');
    
    % wvfGet(wvf,'ref pupil plane size')
    % wvf = wvfSet(wvf,'ref pupil plane size',16.2120);
    wvf = wvfSet(wvf,'ref pupil plane size',50);
    % wvf = wvfSet(wvf,'ref pupil plane size',25);
    
    wvf = wvfComputePSF(wvf);
    
    % This plot is not the same in ISETBio and ISETCam
    % wvfPlot(wvf,'image pupil phase','mm');
    wvfPlot(wvf,'psf space','um',wave,10);
    % wvfPlot(wvf,'image pupil amp','mm')
    
    %% Get the parameters we need for the search
    
    thisWaveUM  = wvfGet(wvf,'wave','um');
    thisWaveNM  = wvfGet(wvf,'wave','nm');
    
    % In ISETCam this is simmply pupil diameter.  In ISETBio this is the
    % calc pupil diameters
    pupilSizeMM = wvfGet(wvf,'calc pupil diameter','mm');
    zpupilDiameterMM = wvfGet(wvf,'measured pupil diameter');
    
    pupilPlaneSizeMM = wvfGet(wvf,'pupil plane size','mm',thisWaveNM);
    nPixels = wvfGet(wvf,'spatial samples');
    wvf     = wvfComputePSF(wvf);
    
    % These are the spatial samples and psf value.
    % When we get a psf from another source, we should interpolate the
    % values to these spatial samples.
    samp      = wvfGet(wvf, 'psf spatial samples', 'um', wave);
    psfTarget = wvfGet(wvf,'psf');   % The spatial sample positions
    
    %{
% This example interpolates the data from t_gullstrandEyeTrace into
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
    disp(x);
    
    % Compare the PSFs to make sure we have a match
    wvf2 = wvfSet(wvf,'zcoeffs',x);
    wvf2 = wvfComputePSF(wvf2);
    wvfPlot(wvf2,'image psf space','um');
    title('Estimated PSF');
    
    vcNewGraphWin;
    imagesc(samp,samp,psfTarget);
    title('Target PSF'); axis image; colormap(hot);
    grid on
    
    % Error in PSF
    psf = wvfGet(wvf2,'psf',wave);
    vcNewGraphWin([],'tall');
    subplot(2,1,1), mesh(psfTarget - psf);
    title('Target - True PSF'); colormap(hot);
    grid on;

    subplot(2,1,2)
    plot(psfTarget(1:5:end),psf(1:5:end),'.');
    axis equal; identityLine;
    xlabel('Target'); ylabel('Estimate');
    
    %%
    fprintf('True zcoeffs\n');
    disp(zcoeffs(1:nCoeffs));
    
    %% Show the pupil phase functions
    %{
    vcNewGraphWin([],'tall');
    subplot(2,1,1)
    wvfPlot(wvf,'image pupil phase','mm',wave,'no window');
    
    wvf2 = wvfSet(wvf,'zcoeffs',x);
    wvf2 = wvfComputePSF(wvf2);
    subplot(2,1,2) 
    wvfPlot(wvf2,'image pupil phase','mm',wave,'no window');
    %}
end

%%