%% Draft for an app to calculate cone fundamentals
%

%% Cone fundamental plots

theCones = coneMosaic;
theCones.wave = 400:1:700;
theLens = Lens;
theLens.wave = 400:1:700;
wave = theCones.wave;

theCones.macular.density = 0.35; % 0.35 is Default
theCones.macular.eccDensity(0:1:15)
theLens.density = 1; % 1 is default

%% 
ieNewGraphWin;
plot(wave,theCones.pigment.quantaFundamentals);

%% Macular pigment and photopigment
ieNewGraphWin;
plot(wave,theCones.qe);

%%
ieNewGraphWin([],'wide');
subplot(1,2,1)
coneFundamentals = bsxfun(@times, theLens.transmittance', theCones.qe);
plot(wave,coneFundamentals,'Linewidth',2);
xlabel('Wavelength (nm)');
ylabel('Efficiency')
grid on;

subplot(1,2,2)
theCones.macular.density = theCones.macular.eccDensity(10); % 0.35 is Default
coneFundamentals = bsxfun(@times, theLens.transmittance', theCones.qe);
plot(wave,coneFundamentals,'Linewidth',2);
xlabel('Wavelength (nm)');
ylabel('Absorptance')
grid on;

%%   Noise distribution graphs

lambda = [1 5 20];
nSamp = 10000;

%% Poisson
ieNewGraphWin;

val = zeros(nSamp,numel(lambda));
for ii=1:numel(lambda)
    val(:,ii) = iePoisson(lambda(ii),'nSamp',[nSamp,1]);
end

h = histogram(val(:,1)); hold on;
h.BinLimits = [-5 40]; h.BinWidth = 1;
h = histogram(val(:,2)); 
h.BinLimits = [-5 40]; h.BinWidth = 1;
h = histogram(val(:,3));
h.BinLimits = [-5 40]; h.BinWidth = 1;
grid on
title('Poisson');

%% Normal

ieNewGraphWin;

val = zeros(nSamp,numel(lambda));
for ii=1:numel(lambda)
    val(:,ii) = sqrt(lambda(ii))*randn(nSamp,1) + lambda(ii);
end

h = histogram(val(:,1)); hold on;
h.BinLimits = [-5 40]; h.BinWidth = 1;
h = histogram(val(:,2)); 
h.BinLimits = [-5 40]; h.BinWidth = 1;
h = histogram(val(:,3));
h.BinLimits = [-5 40]; h.BinWidth = 1;
grid on
title('Normal');

