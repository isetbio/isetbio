function txt = summary(thisEye,varargin)
% Summarize the PBRT scene and rendering parameters
%
% Synopsis
% txt = sceneEye.summary(thisEye,varargin)
% Input
%   thisEye:  sceneEye that includes a recipe
%
% Optional key/val pairs
%
% Output
%   txt:  Text with the information
%
% Description
%
% 

% Examples:
%{
   
%}
%% We will add arguments in the future

thisR = thisEye.get('recipe');
if isempty(thisR), disp('No rendering recipe found'); return; end

%% Basic information printed for now

delimit = '---------------';
txt = addText(delimit,sprintf('\nEye Model: %s\n',thisR.get('camera subtype')));
txt = addText(txt,sprintf('Lens: %s\n',thisR.get('lens file')));
txt = addText(txt,sprintf('Resolution: %d %d\n',thisR.get('spatial resolution')));
txt = addText(txt,sprintf('Rays per pixel: %d \n',thisR.get('rays per pixel')));
txt = addText(txt,sprintf('Ray bounces: %d \n',thisR.get('n bounces')));
txt = addText(txt,sprintf('Scene: %s \n',thisR.get('input basename')));
txt = addText(txt,sprintf('Docker dir: %s \n',thisR.get('output dir')));
txt = addText(txt,delimit);

%% Will get fancier in the future

disp(txt);

end

