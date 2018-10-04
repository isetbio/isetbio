%% t_wvfMTF
%
%  Calculate and plot the MTF for different defocus values on the standard
%  human wvf
% 
%% Show the PSF
wvf = wvfCreate;
wvf = wvfComputePSF(wvf);
wvfPlot(wvf,'psf','min');

%%
wvfD = wvfSet(wvf,'zcoeffs',1,{'defocus'});
wvfD = wvfComputePSF(wvfD);
wvfPlot(wvfD,'psf','min');

%% Read to compute the MTF/OTF now!

%%