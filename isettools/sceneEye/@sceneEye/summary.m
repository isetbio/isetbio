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
if thisEye.get('use pinhole'), renderMode = 'scene'; disp('In pinhole mode'), end

%% Basic information printed for now

delimit = '---------------';
txt = addText(delimit,sprintf('\nHuman Eye Model: %s\n',thisR.get('camera subtype')));
txt = addText(txt,sprintf('Lens: %s\n',thisR.get('lens basename')));
txt = addText(txt,sprintf('Rays per pixel: %d \n',thisR.get('rays per pixel')));
txt = addText(txt,sprintf('Ray bounces: %d \n',thisR.get('n bounces')));
txt = addText(txt,sprintf('Integrator:  %s\n',thisR.get('integrator')));

txt = addText(txt,delimit);
txt = addText(txt,sprintf('\n*** %s *** \n',renderMode));
txt = addText(txt,sprintf('Input: %s \n',thisR.get('input basename')));
txt = addText(txt,sprintf('Docker dir: %s \n',thisR.get('output dir')));

txt = addText(txt,delimit);
if ~thisEye.usePinhole
    % Camera
    focallength = thisR.get('film distance');   % In meters
    txt = addText(txt,sprintf('\nFocal length: %0.4f (m), %0.2f (diopters)',focallength,1/focallength));
    txt = addText(txt,sprintf('\nPupil diameter: %0.1f mm',thisR.get('pupil diameter','mm')));
end

txt = addText(txt,sprintf('\nFOV: %0.1f deg',thisR.get('fov')));
txt = addText(txt,sprintf('\nSpatial samples: %d %d',thisR.get('spatial resolution')));
txt = addText(txt,sprintf('\nLens pigment density: %0.1f',thisEye.get('lens density')));

caEnabled = thisR.get('chromatic aberration');
if caEnabled, nBands = thisR.get('num ca bands');
    txt = addText(txt,sprintf('\nChromatic aberration: Enabled (%d bands)',nBands));
else
    txt = addText(txt,sprintf('\nChromatic aberration: Disabled'));
end

diffEnabled = thisR.get('diffraction');
switch diffEnabled
    case 'true'
        txt = addText(txt,sprintf('\nDiffraction: Enabled\n'));
    case 'false'
        txt = addText(txt,sprintf('\nDiffraction: Disabled\n'));
    otherwise
        txt = addText(txt,sprintf('\nDiffraction: Bad string\n'));
end
txt = addText(txt,delimit);

%% Will get fancier in the future
disp(txt);

end

