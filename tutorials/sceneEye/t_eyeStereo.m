%% t_eyeStereo
% SkipFile
% Stereo sceneEye tutorial uses obsolete constructor/API patterns and
% multiple CPU human-eye renders; keep it interactive for now.
%
% Please read t_eyeIntro first.
%
% Create a stereo pair of retinal irradiance images.
%
% See also
%   t_eyeCrop2Cones
%


%% Make an oi of the chess set scene using the LeGrand eye model

thisSE = sceneEye('chess set scaled','human eye','navarro');

thisSE.set('lens density',0);   % Just because I can

thisSE.set('rays per pixel',512);  % Pretty quick, but not high quality

oiLeft = thisSE.render;  % Render radiance and depth, and then show
oiWindow(oiLeft);

%% Shift the eye position

% Change the eye position (from) but stay focused on the same object (to).
% I shifted the eye position by a lot (12 mm) so the image difference is be
% easy to see.  The inter-pupil difference is really only 6-8 cm 

fromLeft = thisSE.get('from');         % Current camera location
fromRight = fromLeft + [6,0,0]*1e-2;   % Shift it 6 cm
thisSE.set('from',fromRight);  

oiRight = thisSE.render('render type','radiance');
oiWindow(oiRight);

%% END
