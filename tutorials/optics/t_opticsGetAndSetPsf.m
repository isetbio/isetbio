% t_opticsGetAndSetPsf
%
% Show how to get and set psf's in units of minutes of arc, or from otfs
% specified over frequencies in cycles/deg. This ability is useful when we
% want to build isetbio optics using various estimates of optical quality
% in the literature.

% 2/5/17  dhb  Wrote starting with an extant tutorial.

%% Initialize
clear; ieInit;

%% Create the oi structure and pull out the optics
%
% Snag the wavelengths while we're at it.
oi = oiCreate('wvf human');
optics = oiGet(oi,'optics');
wls = opticsGet(optics,'wave');

%% Make the support an even number of samples, if it is odd
%
% The PTB code for optics requires even number of samples, because at one
% point we were more confident of how space and frequency got packed for
% even samples than for odd.  That could be relaxed in the PTB code, but
% the approach we take here is to work in isetbio with support that has an
% even number of samples.
%
% Check is only on x spatial frequencies.  This is OK as long as support is square, 
% which is almost surely true.  We do check that, though.
sfValuesCyclesMm = opticsGet(optics,'otf support','mm');
if (length(sfValuesCyclesMm{1}) ~= length(sfValuesCyclesMm{2}))
    error('Our code assumes that sf support for otf is square, but it isn''t here.')
end
if (rem(length(sfValuesCyclesMm{1}),2) ~= 0)
    % Lop off highest postive frequency in both x and y and put back.
    sfValuesCyclesMm{1} = sfValuesCyclesMm{1}(1:end-1);
    optics = opticsSet(optics,'otffx',sfValuesCyclesMm{1});
    sfValuesCyclesMm{2} = sfValuesCyclesMm{2}(1:end-1);
    optics = opticsSet(optics,'otffy',sfValuesCyclesMm{2});
    
    % Lop off highest postive frequency in the otf. Isetbio stores the otf
    % in Matlab's first entry is zero freuqency format, but the lopping is
    % easier to think about in the zero frequency at the center format.
    % Use fftshift and ifftshift to go back and forth, and do the lopping
    % in between.  We need to do all the wavelengths, and put the whole
    % resized cube back in, in one fell swoop.
    otfOddSupport = opticsGet(optics,'otf data');
    for ii = 1:length(wls)
        otfCentered = fftshift(otfOddSupport(:,:,ii));
        otfCentered = otfCentered(1:end-1,1:end-1);
        otfEvenSupport(:,:,ii) = ifftshift(otfCentered);
    end
    optics = opticsSet(optics,'otf data',otfEvenSupport);
    clear otfOddSupport otfEvenSupport otfCentered
end

%% Get the gridded spatial frequency support of the otf in cycles/deg.
%
% We'll also keep it around in cycles/mm.
%
% And convert to support in cycles per degree using 300 um per degree,
% which is the number that appears to be baked into the optics object.
uMPerMm = 1000;
uMPerDegree = 300;
[xSfGridCyclesMm,ySfGridCyclesMm] = meshgrid(sfValuesCyclesMm{1},sfValuesCyclesMm{2});
xSfGridCyclesDegree = uMPerDegree*xSfGridCyclesMm/uMPerMm;
ySfGridCyclesDegree = uMPerDegree*ySfGridCyclesMm/uMPerMm;

%% Get isetbio format OTF at a specified wavelength
theWl = 550;
otf = opticsGet(optics,'otf data',theWl);

%% Derive the psf from the otf
%
% We have to convert to the zero sf at center representation to use
% OtfToPsf, using ifftshift.
[xGridMinutes,yGridMinutes,psf] = OtfToPsf(xSfGridCyclesDegree,ySfGridCyclesDegree,ifftshift(otf));
centerPosition = length(sfValuesCyclesMm{1})/2+1;
position1DMinutes = xGridMinutes(centerPosition,:);
wvfHuman1DPsf = psf(centerPosition,:);

%% Show otf and psf in 2D
figure; clf;
subplot(1,2,1); hold on;
mesh(xSfGridCyclesDegree,ySfGridCyclesDegree,fftshift(abs(otf)));
xlim([-100 100]); ylim([-100 100]);
axis('square');
xlabel('X SF (cycles/deg)');
ylabel('Y SF (cycles/deg)');
title('OTF');
subplot(1,2,2); hold on
mesh(xGridMinutes,yGridMinutes,psf);
xlim([-10 10]); ylim([-10 10]);
axis('square');
xlabel('x (minutes)');
ylabel('y (minutes)');
title('PSF');

%% Get Davila-Geisler PSF and plot against what isetbio wvf looks like
DavilaGeislerLsf = DavilaGeislerLSFMinutes(position1DMinutes);
DavilaGeislerPsf = LsfToPsf(DavilaGeislerLsf);
psfFig = figure; clf; hold on
plot(position1DMinutes,wvfHuman1DPsf/max(wvfHuman1DPsf),'r','LineWidth',3);
plot(position1DMinutes,DavilaGeislerPsf(centerPosition,:)/max(DavilaGeislerPsf(centerPosition,:)),'g-','LineWidth',2);
xlim([-4 4]);
xlabel('Position (minutes');
ylabel('Normalized PSF Slice');
title('PSF')
legend({sprintf('Wvf Human @%d nm',theWl),'Davila-Geisler'});

%% Stick DavilaGeisler into the optics structure
%
% The ifftshift puts things into the isetbio format.
[~,~,DavilaGeislerOtfCentered] = PsfToOtf(xGridMinutes,yGridMinutes,DavilaGeislerPsf);
DavilaGeislerOtfIsetbio = ifftshift(DavilaGeislerOtfCentered);
insertOtf = zeros(size(opticsGet(optics,'otf data')));
for ii = 1:length(wls)
    insertOtf(:,:,ii) = DavilaGeislerOtfIsetbio;
end
optics = opticsSet(optics,'otf data',insertOtf);

%% Make sure everything is hunky-dory by making the plot using isetbio's fcn\
oi = oiSet(oi,'optics',optics);
udata = oiPlot(oi,'psf',[],theWl);

figure(psfFig);
plot(60*udata.x(centerPosition,:)/300,udata.psf(centerPosition,:)/max(udata.psf(centerPosition,:)),'k:','LineWidth',2);





