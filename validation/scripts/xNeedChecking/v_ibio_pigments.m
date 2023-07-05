%% v_pigments
%
% Show macular and lens pigment curve plots
%
%
% Copyright, ISETBIO Team, 2014

%%
ieInit

%% Plot macular pigment absorptance at a series of densities
% Absorptance is the fraction of photons absorbed

m = Macular;
wave = m.wave;

dList = 0:.1:.5;

vcNewGraphWin;
for ii=1:length(dList)
    m.density = dList(ii);
    sa = m.absorptance;
    plot(wave,sa);
    hold on;
end
xlabel('Wavelength')
ylabel('Spectral absorptance')
title('Macular pigment photon absorptions')
legend(num2str(dList'));

%% Now, plot the spctral density
% This is the unit density scaled by the actual pigment density

vcNewGraphWin;
for ii=1:length(dList)
    m.density = dList(ii);
    sa = m.spectralDensity;
    plot(wave,sa);
    hold on;
end
xlabel('Wavelength')
ylabel('Spectral density')
title('Macular pigment spectral densities')

legend(num2str(dList'));

%%
