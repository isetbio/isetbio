function txt = summary(thisEye,varargin)
% Summarize the PBRT scene and rendering parameters
%
% Synopsis
% txt = sceneEye.summary(thisEye,varargin)
%
% Description
%   Shows the user some of the key parameters used by PBRT and about the
%   scene.  
%
% Input
%   thisEye:  sceneEye that includes a recipe
%
% Optional key/val pairs
%   N/A
%
% Output
%   txt:  Text with the information
%
% See also
%   sceneEye, sceneEye.get

%% We will add arguments in the future

thisR = thisEye.get('recipe');
if isempty(thisR), disp('No rendering recipe found'); return; end

renderMode = 'optical image';
if thisEye.get('use pinhole'), renderMode = 'scene'; end

%% Basic information printed for now

delimit = '---------------';
txt = addText(delimit,sprintf('\nEye Model: %s\n',thisR.get('camera subtype')));
txt = addText(txt,sprintf('Lens: %s\n',thisR.get('lens file')));
txt = addText(txt,sprintf('Rays per pixel: %d \n',thisR.get('rays per pixel')));
txt = addText(txt,sprintf('Ray bounces: %d \n',thisR.get('n bounces')));
txt = addText(txt,sprintf('Integrator:  %s\n',thisR.get('integrator')));

txt = addText(txt,delimit);
txt = addText(txt,sprintf('\nRendering an *** %s *** \n',renderMode));
txt = addText(txt,sprintf('Input: %s \n',thisR.get('input basename')));
txt = addText(txt,sprintf('Docker dir: %s \n',thisR.get('output dir')));
txt = addText(txt,delimit);
txt = addText(txt,sprintf('\nFocal distance: %0.2f (m)\n',thisR.get('focal distance','m')));
caEnabled = thisR.get('chromatic aberration');
if caEnabled, nBands = thisR.get('num ca bands');
    txt = addText(txt,sprintf('Chromatic aberration: Enabled (%d bands)\n',nBands));
else
    txt = addText(txt,sprintf('Chromatic aberration: Disabled\n'));
end
txt = addText(txt,sprintf('FOV: %0.1f deg\n',thisR.get('fov')));
txt = addText(txt,sprintf('Spatial samples: %d %d\n',thisR.get('spatial resolution')));

txt = addText(txt,delimit);


%% Will get fancier in the future

disp(txt);

end

