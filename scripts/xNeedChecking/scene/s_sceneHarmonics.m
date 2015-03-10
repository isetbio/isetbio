% s_sceneHarmonics

%   To create a pattern that is the sum of one or more harmonics,
%   use sceneCreate as in
%   TwoFreq= sceneCreate('harmonic', params);
%   Below, we explain the params data

% Copyright: ImagEval Consulting 2011

%% Creating a scene that has a single harmonic frequency
%   The harmonic call can set the frequency, contrast, phase,
%  angle, row and col size of the harmonic.  
%  The frequency unit in this case is cycles/image.  
%   To obtain cycles per degree, divide by the field of view.

% The parameters that defines the harmonic are:
%        params.freq = 1; params.contrast = 1; params.ph = 0;
%        params.ang= 0; params.row = 128; params.col = 128; 
%        params.GaborFlag=0;

% call sceneCreate with the params
%        [scene,params] = sceneCreate('harmonic',params);   
%

%% To create a pattern that is the sum of one or more harmonics, 

% params.freq, params.contrast, params.ang and params.ph can be vectors 
    
params.freq =  [1 5]; % spatial frequencies of 1 and 5
params.contrast = [0.2, 0.6]; % contrast of the two frequencies
params.ang  = [0, 0]; % orientations
params.ph  = [0 0]; % phase
TwoFreq = sceneCreate('harmonic',params);
vcAddAndSelectObject(TwoFreq);
sceneWindow

%% vary the orientation of the two harmonics
params.freq =  [2 5]; % spatial frequencies of 1 and 5
params.contrast = [0.6, 0.6]; % contrast of the two frequencies
params.ang  = [pi/4, -pi/4]; % orientations
params.ph  = [0 0]; % phase
TwoFreq = sceneCreate('harmonic',params);
vcAddAndSelectObject(TwoFreq);
sceneWindow

%%  another variation
params.freq =  [5 5]; % spatial frequencies of 1 and 5
params.contrast = [0.6, 0.6]; % contrast of the two frequencies
params.ang  = [pi/4, -pi/4]; % orientations
params.ph  = [0 0]; % phase
TwoFreq = sceneCreate('harmonic',params);
vcAddAndSelectObject(TwoFreq);
sceneWindow