% t_oisCreate
%
% Create oiSequences using the oisCreate() function.  This functions
% provides a general interface for a few classic psychophysical stimuli.
%
% The oiSequence() class is illustrated in t_oiSequence.m
%
% BW, ISETBIO TEAM, 2017

%%
ieInit

%% Harmonic
%
% This code relies on the sceneCreate('harmonic', ...) function to create
% harmonics.  Typically, we create a pair where one has positive contrast
% and a second has 0 contrast.  This produces a simple modulating oi
% sequence

clear hparams
% The modulating harmonic parameters.  The possibilities are defined in the
% sceneCreate('harmonic', ...) function, and the complete list of
% parameters is returned by harmonicP.
hparams(2) = harmonicP; 
hparams(2).freq = 6; hparams(2).GaborFlag = .2;

% The matched, zero contrast, harmonic parameters
hparams(1) = hparams(2); hparams(1).contrast = 0;

% The general scene properties can also be set.  Here we make the scene 1
% deg of visual angle
sparams.fov = 1;

% And then we make a Gaussian temporal modulation that brings the stimulus
% on and off
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);

% The two harmonics are 'blended', which means at each moment in time we
% have a weighted sum of the two where the weights sum to 1.
ois = oisCreate('harmonic','blend',stimWeights, 'testParameters',hparams,'sceneParameters',sparams);

ois.visualize;

%% We can change the parameters for different effects

% Frequency changin
hparams(1) = hparams(2); hparams(1).freq = 4;
ois = oisCreate('harmonic','blend',stimWeights, 'testParameters',hparams,'sceneParameters',sparams);
ois.visualize;

%% A little phase shifting

hparams(1) = hparams(2); hparams(1).ph = 0;
ois = oisCreate('harmonic','blend',stimWeights, 'testParameters',hparams,'sceneParameters',sparams);
ois.visualize;

%% Vernier
%
% The line offset case is managed using the sceeCreate('vernier', ...)
% scene generations.  The code for the vernier stimuli needs
% simplification, checking, and improvement. (BW).
% 

clear vparams; 
% The vernier parameter defaultss
vparams(2) = vernierP;

% A black background
vparams(2).name = 'offset'; vparams(2).bgColor = 0; 

% The bar is added to the black background
vparams(1) = vparams(2);
vparams(1).barWidth = 0; vparams(1).bgColor = 0.5; vparams(1).name = 'uniform';

% The scene
sparams.fov = 1;

% The temporal modulation
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);

% You can return the scenes, if you like.
[vernier, scenes] = oisCreate('vernier','add', stimWeights,'testParameters',vparams,'sceneParameters',sparams);
vernier.visualize;

% To see the scenes, prior to creating the optical image
% ieAddObject(scenes{1}); ieAddObject(scenes{2}); sceneWindow;

%% Change up the parameters
vparams(2).barColor = 0.5; 

vernier = oisCreate('vernier','add', stimWeights,'testParameters',vparams,'sceneParameters',sparams);
vernier.visualize;



%% Impulse (temporal)

clear iparams

sparams.fov = 1; sparams.luminance = 100;
stimWeights = zeros(1,50); stimWeights(10:13) = 1;

impulse = oisCreate('impulse','add', stimWeights,'sceneParameters',sparams);

impulse.visualize;

%%