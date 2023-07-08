%% v_ibio_wvfZernikePolynomials
%
% Make plots of various pupil functions and their respective point-spread
% functions for different Zernike polynomials of 1st through 3rd radial
% orders (OSA j indices 1 through 9).
%
% Each time through the loop we see the effect of wiggling one coefficient.
%
% In this loop, we plot the wavefront aberrations (measured in microns),
% the pupil function phase (radians), and the PSF.
%
% The wavefront aberration plots we get match those
%  http://www.telescope-optics.net/monochromatic_eye_aberrations.htm
% except for defocus, where we have opposite sign.

%%
wvf0 = wvfCreate;
wvf0 = wvfSet(wvf0,'calculated pupil',wvfGet(wvf0,'measured pupil','mm'));
pupilfuncrangeMM = 4;
wList = wvfGet(wvf0,'calc wave');
jindices = 0:9;
maxMM = 4; 
for ii = jindices
    vcNewGraphWin([],'tall');
    insertCoeff = 1;
    wvf = wvfSet(wvf0,'zcoeffs',insertCoeff,ii);
    wvf = wvfCompute(wvf);
    [n,m] = wvfOSAIndexToZernikeNM(ii);

    subplot(3,1,1);
    wvfPlot(wvf,'2d wavefront aberrations space','mm',[],pupilfuncrangeMM,'no window');
    title(sprintf('Wavefront aberrations for j = %d (n = %d, m = %d)',ii,n,m));

    subplot(3,1,2);
    wvfPlot(wvf,'2d pupil phase space','mm',wList,pupilfuncrangeMM,'no window');
    title(sprintf('Pupil function phase for j = %d (n = %d, m = %d)',ii,n,m));

    subplot(3,1,3);
    wvfPlot(wvf,'2d psf space','mm',wList,maxMM,'no window');
    
    % Save out what it does today
    %
    % Little bit of nonsense here to avoid a "-" in the identifier string,
    % yet still distinguish plus from minus.
    if (m < 0)
        mStr = ['m' num2str(abs(m))];
    else
        mStr = num2str(m);
    end
    if (n < 0)
        nStr = ['m' num2str(abs(n))];
    else
        nStr = num2str(n);
    end
    
end
