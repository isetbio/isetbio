% s_oisHarmonic
%
% Create an oiSequence of a harmonic
%

%% 
ieInit
clear params
clear scene
%%

scene = cell(1,2);

params.freq =  10; % spatial frequencies of 1 and 5
params.contrast = 0.6; % contrast of the two frequencies
params.ang  = [0, 0]; % orientations
params.ph  = [0 0]; % phase
scene{1} = sceneCreate('harmonic',params);
scene{1} = sceneSet(scene{1},'name',sprintf('F %d',params.freq));
ieAddObject(scene{1});

% Create scene: background field only, no harmonic
clear params
params.freq =  0; % spatial frequencies of 1 and 5
params.contrast = 0; % contrast of the two frequencies
params.ang  = [0, 0]; % orientations
params.ph  = [0 0]; % phase
scene{2} = sceneCreate('harmonic',params);
scene{2} = sceneSet(scene{2},'name','Uniform');
ieAddObject(scene{2});

%%
imgFov = .5 ;      % image field of view
vDist  = 0.3;          % viewing distance (meter)

% set scene fov
for ii = 1 : length(scene)
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end
% sceneWindow

%%

% These are the default.
oi = oiCreate('wvf human');

% Compute optical images from the scene
OIs = cell(length(scene), 1);
for ii = 1 : length(OIs)
    OIs{ii} = oiCompute(oi,scene{ii});
end


%% Build the oiSequence

% We build the stimulus using a time series of weights. We have the mean
% field on for a while, then rise/fall, then mean field.
zTime = 50;   % Mean field beginning and end (ms)
stimWeights = fspecial('gaussian',[1,50],15);
stimWeights = ieScale(stimWeights,0,1);
weights = [zeros(1, zTime), stimWeights, zeros(1, zTime)];

% Temporal samples.  Typically 1 ms, which is set by the parameter in the
% cone mosasic integration time.  That time is locked to the eye movements.
tSamples = length(weights); 
sampleTimes = 0.002*(1:tSamples);  % Time in sec

% vcNewGraphWin; plot(1:tSamples, weights,'o');
% xlabel('Time (ms)');

% The weights define some amount of the constant background and some amount
% of the line on the same constant background
oiHarmonicSeq = oiSequence(OIs{2}, OIs{1}, ...
    sampleTimes, weights, ...
    'composition', 'blend');
oiHarmonicSeq.visualize('format','movie');

%%