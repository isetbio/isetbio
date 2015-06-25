%% s_wvf2OIHumanSamples
%
% Check the variation in the Thibos wavefront model.
%
% See also: s_wvf2OIHuman and s_wvf2OI
%
% 7/23/12  dhb  Drop off what I believe to be j=0 value from Thibos model.
%
% (c) Wavefront Toolbox Team, 2012


%% Initialize
s_initISET

%%
maxUM = 10;
wave = 400:10:700; wave = wave(:);
pupilMM = 3;

scene = sceneCreate; vcAddObject(scene);
oi = oiCreate('human'); vcAddObject(oi);

%%  Create some examples.
% Can use either Thibos statistical model, or read in the measurements we
% got from Heidi Hofer
whichTypeOfSamples = 'ThibosStatiscalModel';
switch (whichTypeOfSamples)
    case 'ThibosStatiscalModel'
        N = 10;
        [sample_mean S] = wvfLoadThibosVirtualEyes(pupilMM);
        zSamples = ieMvnrnd(sample_mean,S,N)';
        measPupilSizeMM = pupilMM;
        measWavelengthNM = 550;
        
    case 'HoferMeasurements'
        zSamples = importdata('sampleZernikeCoeffs.txt');
        if (size(zSamples,1) ~= 66 || size(zSamples,2) ~= 9)
            error('Surprising size for read in Hofer zSamples.')
        end
        N = size(zSamples,2);
        measPupilSizeMM = 6;
        measWavelengthNM = 550;
        
    otherwise
        error('Unknown data source specified');
end
nCoeffs = size(zSamples,1);


%% Convert WVF human data to ISET
oiD = cell(N,1);

% Create samples
for ii=1:N
    name = sprintf('%d human-%d',ii,pupilMM);
    wvfP = wvfCreate('wave',wave,'name',name);
    wvfP = wvfSet(wvfP,'measured wavelength',measWavelengthNM);
    wvfP = wvfSet(wvfP,'measured pupil size',measPupilSizeMM);
    wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
    
    z = wvfGet(wvfP,'zcoeffs');
    z(1:nCoeffs) = zSamples(:,ii);
    wvfP = wvfSet(wvfP,'zcoeffs',z);
    
    wvfP = wvfComputePSF(wvfP);
    oiD{ii} = wvf2oi(wvfP,'human');
    oiD{ii} = oiSet(oiD{ii},'name',name);
end

%% Now compare the slanted bar response in the OI
% These are reasonably close for calculations separated by so many years.

scene  = sceneCreate('slanted bar');
scene  = sceneSet(scene,'h fov',1);
bb = blackbody(sceneGet(scene,'wave'),6500,'energy');
scene = sceneAdjustIlluminant(scene,bb);
% vcAddAndSelectObject(scene); sceneWindow;

sensor = sensorCreate('human');
sensor = sensorSet(sensor,'noise flag',0);   % Turn off all noise.
sensor = sensorSet(sensor,'exp time',0.050);
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'hfov'),scene,oiD{1});

%%
uData = cell(N,1);
for ii=1:N
    oiD{ii} = oiCompute(oiD{ii},scene);
    sensorD = sensorCompute(sensor,oiD{ii});
    sensorD = sensorSet(sensorD,'name','Thibos calc');
    vcAddAndSelectObject(sensorD); sensorWindow('scale',1);
    [uData{ii},g] = sensorPlot(sensorD,'electrons hline',[1,80]);
    close(g)
end

%%
vcNewGraphWin([],'tall');
hold on
c = {'r-','g-','b-'};
for ii=1:N
    for jj=1:3
        hold on; grid on;
        subplot(3,1,jj)
        plot(uData{ii}.pos{jj},uData{ii}.data{jj},c{jj})
        title(sprintf('%s, pupil size %0.1f mm',whichTypeOfSamples,pupilMM));
    end
end

%%
% sensorMW = sensorCompute(sensor,oiMW);
% sensorMW = sensorSet(sensorMW,'name','MW calc');
% vcAddAndSelectObject(sensorMW); sensorWindow('scale',1);
% sensorPlot(sensorMW,'electrons hline',[1,80]);
% title('Marimont and Wandell')

%% End
