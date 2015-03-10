%s_opticsDLPSF
%
% In progress - not running properly yet.
%
% Plots of the diffraction-limited (DL) PSF as a function of position and
% angle.
%
% A main point here is that the spread as a function of angle is
% independent of f-number
%
% The diffraction limited spread grows in space.
%
% This script shows two ways of calculating the angle in radians.
%
% (c) Imageval Consulting, LLC, 2012

%%
ieInit
close all;  % I think this should be in ieInit

%% Create optical image
oi     = oiCreate;

fLength = 0.017;   % Meters
fNumber = 17/3;    % Dimensionless
oi    = oiSet(oi,'optics flength',fLength);
oi    = oiSet(oi,'optics fnumber',fNumber);

%%  Show linespread at different wavelengths
uData   = plotOI(oi,'ls wavelength');
fNumber = oiGet(oi,'optics fnumber','mm');
fLength = oiGet(oi,'optics focal length','mm');
title(sprintf('F/# = %.2f  Foc Leng = %.2f (mm)',fNumber,fLength));

%% Show the linespread in units of arc min rather than position (um)
vcNewGraphWin
posMM    = uData.x/1000;                             % Microns to mm
aMinutes = rad2deg(atan2(posMM,fLength),'arcmin');   % Angle in radians

mesh(aMinutes,uData.wavelength,uData.lsWave);
view(30,20);
xlabel('angle (arc min)');
ylabel('wavelength (nm)');

%% For a wavelength, show the full psf near the Airy Ring
thisWave = 400;
uData = plotOI(oi,'psf',[],thisWave); 

view(2)
AiryRingUM = (2.44*(thisWave/1000)*fNumber);
set(gca,'xlim',[-AiryRingUM AiryRingUM],'ylim',[-AiryRingUM AiryRingUM])

%% Show a slice through the psf as a function of angle
r = size(uData.x,1);
mid = ceil(r/2);
psfMid = uData.psf(mid,:);
posMM  = uData.x(mid,:)/1000;               % Microns to mm
posMinutes = rad2deg(atan2(posMM,fLength),'arcmin');

vcNewGraphWin;
plot(posMinutes,psfMid)
xlabel('Arc min')
AiryRingMM = AiryRingUM/1000;
AiryRingMinutes = rad2deg(atan2(AiryRingMM,fLength),'arcmin'); % Radians
set(gca,'xlim',2*[-AiryRingMinutes AiryRingMinutes])

pDiameter = oiGet(oi,'optics pupil diameter','mm');
str = sprintf('PSF cross section: Wave %d (nm), fNumber %.2f Pupil %.1f (mm)',...
    thisWave,fNumber,pDiameter);
title(str);

%% Or, here is a single line spread
uData = plotOI(oi,'lswavelength');
posMM = uData.x/1000;
aRadians = atan2(posMM,fLength);    % This is angle in radians
aMinutes = rad2deg(aRadians,'arcmin');          % This is angle in arc min
plot(aMinutes,uData.lsWave(1,:),'-',...
    aMinutes,uData.lsWave(16,:),'r:',...
    aMinutes,uData.lsWave(31,:),'g--')
legend('400nm','550nm','700nm')
title('Line spread')
grid on
xlabel('arc min')

%% For wvf comparison, here is just at 550nm
thisWave = 550;
uData = plotOI(oi,'psf',[],thisWave); 

view(2)
AiryRingUM = (2.44*(thisWave/1000)*fNumber);
set(gca,'xlim',[-AiryRingUM AiryRingUM],'ylim',[-AiryRingUM AiryRingUM])

r = size(uData.x,1);
mid = ceil(r/2);
psfMid = uData.psf(mid,:);
posMM = uData.x(mid,:)/1000;               % Microns to mm
posMinutes = rad2deg(atan2(posMM,fLength),'arcmin');

vcNewGraphWin
plot(posMinutes,psfMid)
xlabel('Arc min')
set(gca,'xlim',[-2 2])

%% Now change the f-number and plot in angle again
fNumber = 5*(17/3);
oi = oiSet(oi,'optics fnumber',fNumber);

%%  Show linespread at different wavelengths
plotOI(oi,'ls wavelength');
fNumber = oiGet(oi,'optics fnumber','mm');
fLength = oiGet(oi,'optics focal length','mm');
title(sprintf('F/# = %.2f  Foc Leng = %.2f (mm)',fNumber,fLength));

%%
fNumber = 0.5*(17/3);
oi = oiSet(oi,'optics fnumber',fNumber);

%%  Show linespread at different wavelengths
uData = plotOI(oi,'ls wavelength');
fNumber = opticsGet(optics,'fnumber','mm');
fLength = opticsGet(optics,'focal length','mm');
title(sprintf('F/# = %.2f  Foc Leng = %.2f (mm)',fNumber,fLength));

%% The psf spread in distance only cares about f number
% We can change the f length and the functions is the same
oi = oiSet(oi,'optics flength',fLength);
uData1 = plotOI(oi,'psf 550');

optics = opticsSet(optics,'flength',fLength*10);
oi = oiSet(oi,'optics',optics);
uData2 = plotOI(oi,'psf 550');

%% For a fixed aperture, the spread in angle is invariant to f length.
% So, if we change the fnumber, and plot w.r.t angle, we get the same
% spread function.

% TEST THIS HERE.


aphi = sqrt(aRadians(:,1).^2 + aRadians(:,2).^2);
vcNewGraphWin([],'tall');
subplot(3,1,1), imagesc(aphi); axis image; colorbar
subplot(3,1,2), plot(phi10(:),aphi(:),'.'); axis square; grid on
subplot(3,1,3), plot(phi(:),aphi(:),'.');axis square; grid on

%% Extra

sSupport = oiGet(oi,'spatial support','mm');  % row,col, x/y
dist = sqrt(sSupport(:,:,1).^2 + sSupport(:,:,2).^2);
optics = oiGet(oi,'optics');
fLength  = opticsGet(optics,'focal length','mm');

% tan(phi) = Opp/adjacent; phi = atan2(opp,adjacent);
phi = atan2(dist,fLength);
vcNewGraphWin([],'tall');
subplot(2,1,1), imagesc(phi); axis image; colorbar
subplot(2,1,2), imagesc(dist); axis image; colorbar
fNumber = opticsGet(optics,'f number');

%%
optics = opticsSet(optics,'fnumber',fNumber*2);
oi = oiSet(oi,'optics',optics);
sSupport = oiGet(oi,'spatial support','mm');  % row,col, x/y
dist = sqrt(sSupport(:,:,1).^2 + sSupport(:,:,2).^2);
fLength  = opticsGet(optics,'focal length','mm');
% tan(phi) = Opp/adjacent; phi = atan2(opp,adjacent);
phi2 = atan2(dist,fLength);
vcNewGraphWin([],'tall'); 
subplot(2,1,1), imagesc(phi2); axis image; colorbar
subplot(2,1,2), plot(phi(:),phi2(:),'.');

%%
optics = opticsSet(optics,'fnumber',fNumber*10);
oi = oiSet(oi,'optics',optics);
sSupport = oiGet(oi,'spatial support','mm');  % row,col, x/y
dist = sqrt(sSupport(:,:,1).^2 + sSupport(:,:,2).^2);
fLength  = opticsGet(optics,'focal length','mm');
% tan(phi) = Opp/adjacent; phi = atan2(opp,adjacent);
phi10 = atan2(dist,fLength);
vcNewGraphWin([],'tall'); 
subplot(2,1,1), imagesc(phi); axis image; colorbar
subplot(2,1,2), plot(phi(:),phi10(:),'.');
