%% Draft for an app to calculate cone fundamentals
%

%%
theCones = coneMosaic;
theLens = Lens;
wave = theCones.wave;

theCones.macular.density = 0.4 % 0.35 is Default
theCones.macular.eccDensity(0:1:15)
theLens.density = 1; % 1 is default

%% 
ieNewGraphWin;
plot(wave,theCones.pigment.quantaFundamentals);

%% Macular pigment and photopigment
ieNewGraphWin;
plot(wave,theCones.qe);

%%
ieNewGraphWin;
coneFundamentals = bsxfun(@times, theLens.transmittance, theCones.qe);
plot(wave,coneFundamentals);
