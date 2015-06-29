%% p_Toronto_optics
%

%% Optics section
talkD = fullfile(isetbioRootPath,'scripts','toronto');

%% Diffraction limited simulation properties
oi = oiCreate;
oiPlot(oi,'otf',[],550);
oiPlot(oi,'otf',[],450);

%% Effect of fnumber
optics = opticsCreate('one inch');

optics = opticsSet(optics,'fnumber',2);
oi = oiSet(oi,'optics',optics);
oiPlot(oi,'lsf wavelength');
set(gca,'xlim',[-10 10],'xtick',-10:2:10)
title(sprintf('Pupil radius %.1f (mm)',opticsGet(optics,'pupil radius','mm')))

optics = opticsSet(optics,'fnumber',8);
oi = oiSet(oi,'optics',optics);
oiPlot(oi,'lsf wavelength');
set(gca,'xlim',[-10 10],'xtick',-10:2:10)
opticsGet(optics,'pupil radius','mm')
title(sprintf('Pupil radius %.1f (mm)',opticsGet(optics,'pupil radius','mm')))

%%
oi = oiCreate('human');
oiPlot(oi,'lsf wavelength');
oiPlot(oi,'otf',[],600);
oiPlot(oi,'otf',[],420);
% Shift-invariant representation

% Human chromatic aberration example
% Linespread functions by wavelength

% Depth of field example with Maya/Rendertoolbox example

%% Human Line Spread examples (after s_humanLSF.m)

imSize = 128;
scene = sceneCreate('lined65',imSize);     % D65 SPD for a thin line
scene = sceneSet(scene,'fov',0.4);  % Small field of view 
vcReplaceObject(scene);
% sceneWindow;

oi = oiCreate;
oi = oiSet(oi,'spectrum',sceneGet(scene,'spectrum'));
optics = opticsCreate('human');  % Set up for human optics
oi = oiSet(oi,'optics',optics);
oi = oiCompute(scene,oi);
vcReplaceObject(oi);
oiWindow;

%% Make OIs for several different wavelengths

% grid spacing on image (microns)
gSpacing = 20;

scene410 = sceneInterpolateW(scene,410);
oi = oiCompute(scene410,oi);
oi = oiSet(oi,'name','line-410');
vcAddAndSelectObject(oi);
oiPlot(oi,'irradiance image with grid',[],gSpacing)

scene550 = sceneInterpolateW(scene,550);
oi = oiCompute(scene550,oi);
oi = oiSet(oi,'name','line-550');
vcAddAndSelectObject(oi);
oiPlot(oi,'irradiance image with grid',[],gSpacing)

scene690 = sceneInterpolateW(scene,690);
oi = oiCompute(scene690,oi);
oi = oiSet(oi,'name','line-690');
vcAddAndSelectObject(oi);
oiPlot(oi,'irradiance image with grid',[],gSpacing)

oiWindow;

%%  Illustrate the yellow-blue bar in the cone capture level
% rgbFile = 'C:\Users\Wandell\Documents\Talks\20110715 Toronto\images\yellow-blue-LMS\yellow-blue-comparison.png';
rgbFile = fullfile(talkD,'yellow-blue-comparison.png');
scene = sceneFromFile(rgbFile,'rgb',100,'lcdExample.mat');
scene = sceneSet(scene,'fov',1);
scene = sceneAdjustIlluminant(scene,'equalPhotons.mat');
vcAddAndSelectObject(scene); sceneWindow; % display sceneWindow

oi = oiCreate('human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi); oiWindow; % display sceneWindow

sensor = sensorCreate('human');
sensor = sensorSet(sensor,'expTime',0.200);
sensor = sensorSet(sensor,'size',[192,256]);
sensor = sensorCompute(sensor,oi);

vcAddAndSelectObject(sensor); sensorImageWindow; % display sceneWindow


%%  Show plots as the LMS contrasts as the bar spatial frequency changes
%

% There are 3 cycles/per image so the bar squarewave frequencies are about
% 3/fov cpd.
fov = [0.5,1,2,3];
saveFlag = 1;
for ii=1:length(fov)
    scene = sceneSet(scene,'fov',fov(ii));
    oi = oiCompute(scene,oi);
    sensor = sensorCompute(sensor,oi);
    xy = [1, round(sensorGet(sensor,'rows')*0.8)];
    
    f = sensorPlotLine(sensor,'h','electrons','space',xy);

    if saveFlag
        fname = fullfile(talkD,'images','yellow-blue-LMS',sprintf('yellowBlue-%.1f.fig',fov(ii)));
        set(f,'Position',[1250 60 560 932]);
        saveas(f,fname,'fig');
    end
end


%%  Show harmonic contrast LMS plots as the spatial frequency changes
%
oi = oiCreate('human');

sensor = sensorCreate('human');
sensor = sensorSet(sensor,'expTime',0.200);
sensor = sensorSet(sensor,'size',[192,256]);

parms.contrast = 1; parms.ph = 0;
parms.ang= 0; parms.row = 128; parms.col = 128; 
parms.GaborFlag=0;

%% Run the loop Cone modulations

% Adjust mean luminance and frequency
lum = [100, 50, 30, 10];       % cd/m2.
freq = [1 4 8 16 32 64];   %Cycles per deg or per scene?
saveFlag = 0;              % Save to file
for jj=1:length(lum)
    thisLum = lum(jj);
    for ii=1:length(freq)
        parms.freq = freq(ii);
        [scene,parms] = sceneCreate('harmonic',parms);
        scene = sceneSet(scene,'fov',1);  % freq is in cpd now
        scene = sceneAdjustLuminance(scene,thisLum);
        % vcAddAndSelectObject(scene); sceneWindow; 
        oi = oiCompute(scene,oi);
        sensor = sensorCompute(sensor,oi);
        xy = [1, round(sensorGet(sensor,'rows')*0.8)];
        
        [f, uData] = sensorPlotLine(sensor,'h','electrons','space',xy);
        % mx = 0;
        % for aa=1:3, mx = max(mx,max(uData.data{aa}(:))); end
        % mx = 1000*ceil(mx/1000);
        
        axisList = findobj(f,'type','axes');
        for aa=1:length(axisList),
            set(axisList(aa),'ylim',[0 3000],'xlim',[-120 120]);
        end
        
        if saveFlag
            fname = fullfile(talkD,'images','CSF',sprintf('CSF-%.0f-%.0f.fig',freq(ii),thisLum));
            set(f,'Position',[1250 60 560 932]);
            saveas(f,fname,'fig');
        end

    end
end

%% Assemble the L and M cones at the highest frequency into a single curve

lum = 100;
freq = [8 16 32];     %Cycles per deg or per scene?
saveFlag = 0;       % Save to file
uData = cell(length(lum),length(freq));

for jj=1:length(lum)
    thisLum = lum(jj);
    for ii=1:length(freq)
        parms.freq = freq(ii);
        [scene,parms] = sceneCreate('harmonic',parms);
        scene = sceneSet(scene,'fov',1);  % freq is in cpd now
        scene = sceneAdjustLuminance(scene,thisLum);

        oi = oiCompute(scene,oi);
        sensor = sensorCompute(sensor,oi);
        xy = [1, round(sensorGet(sensor,'rows')*0.8)];
        
        [f, tmp] = sensorPlotLine(sensor,'h','electrons','space',xy);
        axisList = findobj(f,'type','axes');
        for aa=1:length(axisList),
            set(axisList(aa),'ylim',[0 3000],'xlim',[-120 120]);
        end
        % Unscaled
        % vcNewGraphWin;
        % plot(tmp.pos{1},tmp.data{1},'ro'); hold on
        % plot(tmp.pos{2},tmp.data{2},'go'); hold on
        % plot(tmp.pos{3},tmp.data{3},'bo'); hold on

        % Scale the level of the different cone classes to modulate around 1
        for cc = 1:2, tmp.data{cc} = tmp.data{cc}/mean(tmp.data{cc}(:)); end
        pos = cat(1,tmp.pos{1}, tmp.pos{2});
        val = cat(1,tmp.data{1}, tmp.data{2});     
        
        [newPos,idx] = sort(pos); newVal = val(idx);
        vcNewGraphWin; plot(newPos,newVal,'-')
        set(gca,'xlim',[-120 120]);
        name = sprintf('freq %.0f',freq(ii));
        set(gcf,'name',name);
        newValFFT = fftshift(abs(fft(newVal - mean(newVal(:)))));
        vcNewGraphWin; plot(newValFFT); set(gcf,'name',name);
        
    end
end

%% Do a depth of field calculation for the goblet scene

