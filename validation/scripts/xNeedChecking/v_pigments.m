%% v_pigments
%
% Show macular and lens pigment curve plots
%
%
% Copyright, ISETBIO Team, 2014

%%
s_initISET

%% Plot macular pigment absorptance at a series of densities
% This is the fraction of absorbed photons

m = macularCreate;
wave = macularGet(m,'wave');

dList = 0:.1:.5;

vcNewGraphWin;
for ii=1:length(dList)
    m = macularSet(m,'density',dList(ii));
    sa = macularGet(m,'absorptance');
    plot(wave,sa);
    hold on;
end
xlabel('Wavelength')
ylabel('Spectral absorptance')
title('Macular pigment photon absorptions')
legend(num2str(dList'));

%% Now, plot the densities
vcNewGraphWin;
for ii=1:length(dList)
    m = macularSet(m,'density',dList(ii));
    sa = macularGet(m,'spectraldensity');
    plot(wave,sa);
    hold on;
end
xlabel('Wavelength')
ylabel('Spectral absorbance')
title('Macular pigment spectral densities')

legend(num2str(dList'));

%%
