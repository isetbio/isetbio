%% s_sceneIlluminant
%
% Illustrate and test illuminant structure creation and plotting.
%
% (c) Imageval Consulting, LLC 2012

%% Create a blackbody illuminant structure 5000 deg Kelvin

illum = illuminantCreate('blackbody');
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%% Create a blackbody illuminant structure 3000 deg Kelvin
wave = 400:1:700;
illum = illuminantCreate('blackbody',wave,3000);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%% D65, 200 cd/m2
luminance = 200;
wave = [];  % Defaults
illum = illuminantCreate('d65',wave,luminance);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%%
illum = illuminantCreate('equal energy',[],luminance);
e     = illuminantGet(illum,'energy');
w     = illuminantGet(illum,'wave');
mn    = mean(e(:));

vcNewGraphWin; 
plot(w,e); grid on
xlabel('Wavelength'); ylabel('Energy')
set(gca,'ylim',[0.9*mn,1.1*mn]);

%%
illum = illuminantCreate('equal photons',[],luminance);
p     = illuminantGet(illum,'photons');
e     = illuminantGet(illum,'energy');

w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('photons')

vcNewGraphWin; 
plot(w,e); grid on
xlabel('Wavelength'); ylabel('Energy')
%%
illum = illuminantCreate('illuminant C',[],luminance);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')
%%
illum = illuminantCreate('555 nm',[],luminance);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%%
illum = illuminantCreate('d65',400:2:600,200);
e     = illuminantGet(illum,'energy');
w     = illuminantGet(illum,'wave');

vcNewGraphWin([],'tall'); 
plot(w,e,'ks-'); grid on
xlabel('Wavelength'); ylabel('Energy')

%% Now interpolate and overlay.
illum2 = illuminantSet(illum,'wave',400:5:700);
e     = illuminantGet(illum2,'energy');
w     = illuminantGet(illum2,'wave');

hold on;
plot(w,e,'ro'); grid on
xlabel('Wavelength'); ylabel('Energy')

%%
illum = illuminantCreate('fluorescent',400:5:700,10);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%%
illum = illuminantCreate('tungsten',[],300);
p     = illuminantGet(illum,'photons');
w     = illuminantGet(illum,'wave');

vcNewGraphWin; 
plot(w,p); grid on
xlabel('Wavelength'); ylabel('Photons')

%% End