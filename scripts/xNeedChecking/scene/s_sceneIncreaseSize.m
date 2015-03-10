%% v_sceneIncreaseSize
%
%  Illustrate how to increase the scene image size
%
% (c) Imageval Consulting, LLC, 2012

%%
s_initISET

%%
scene = sceneCreate;
vcAddAndSelectObject(scene); sceneWindow;

%% Double the rows, triple the columns
s = [2,3];
p = sceneGet(scene,'photons');
p = imageIncreaseImageRGBSize(p,s);
scene = sceneSet(scene,'photons',p);
vcAddAndSelectObject(scene); sceneWindow;

%% Double the cols
s = [1,2];
p = sceneGet(scene,'photons');
p = imageIncreaseImageRGBSize(p,s);
scene = sceneSet(scene,'photons',p);
vcAddAndSelectObject(scene); sceneWindow;

%% Triple the cols - should return to original aspect ratio

s = [3,1];
p = sceneGet(scene,'photons');
p = imageIncreaseImageRGBSize(p,s);
scene = sceneSet(scene,'photons',p);
vcAddAndSelectObject(scene); sceneWindow;

%%