%% t_imageMultiview
%
% Illustrate how to display multiple RGB images from one of the current
% scene/oi/sensor/vci. 
%
% Copyright Imageval LLC, 2013

%% Initialize
ieInit;

%% To start debugging I ran s_imageIlluminantCorrection
%
% This provides windows with multiple examples
% It takes a little while.
s_sceneChangeIlluminant
%%  Get a list of the objects

% This example shows several of the scene images.  No user interaction
% required.
objType = 'scene';
whichObj = [3 4 5];
% Show in a single window
imageMultiview(objType,whichObj,true);

%% This one allows you to select which ones you want to compare
% They come up in separate windows.
imageMultiview('scene');

%%  To run it in another window, when you have images there, use ...
% sensorWindow;
% imageMultiview('sensor',[ 1 4]);
% imageMultiview('vci',[1 4]);

%% End
