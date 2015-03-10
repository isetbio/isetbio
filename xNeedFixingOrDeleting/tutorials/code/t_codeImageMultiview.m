%% t_imageMultiview
%
% Illustrate how to display multiple RGB images from the scene/oi/sensor/vci objects.
%
% We fill up the various GUI windows with several examples.  Then we show
% how to bring up the RGB images into separate windows.
%
% NOTES:
%  1) Broken because it runs a script that does not currently exist.
%
% Copyright Imageval LLC, 2013

%% Initialize
ieInit;

%% To start debugging I ran s_imageIlluminantCorrection
%
% This provides  windows with multiple examples
% It takes a little while.
s_imageIlluminantCorrection

%%  Get a list of the objects

% This example shows several of the scene images.  No user interaction
% required.
objType = 'scene';
whichObj = [2 3 5];
imageMultiview(objType,whichObj);

%% This one allows you to select which ones you want to compare
oiWindow;
imageMultiview('oi');

%%
sensorWindow;
imageMultiview('sensor',[ 1 4]);

%%
vcimageWindow;
imageMultiview('vci',[1 4]);

%% End
