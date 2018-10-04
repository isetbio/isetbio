%% t_wvfMTF
%
%  Calculate and plot the MTF for different defocus values on the standard
%  human wvf
% 
%%
wvf = wvfCreate;
wvf = wvfComputePSF(wvf);
wvfPlot(wvf,'psf');

