%% t_coneWindowAppearance
%
% Controlling window properties in Matlab.
%
% It is possible to change some of the GUI appearance properties from the
% command line.  Here we illustrate how to make windows invisible, change
% the white point, and change the font size.
%
% NOTES:
%   1) This might be elaborated a bit, to show explicitly more things that
%   you might want to do.
%
% Copyright ImagEval Consultants, LLC, 2010

%% Initialize
ieInit;

%% Opening and closing windows from the command line

% Create an isetbio window
scene = sceneCreate; vcAddAndSelectObject(scene);

sceneWindow('visible','on')
drawnow
pause(1)

sceneWindow('visible','off')
drawnow
pause(1)

%% Interacting with the scene handles
%
% It is also possible to make adjustments to the display by interacting
% with the Matlab handle graphics.  To get the handle to the scene figure,
% you can run
sceneF = ieSessionGet('scene figure');

% The variable sceneF is the handle to the figure
get(sceneF)

% The guidata are available here
guidata(sceneF)

% Or you can get the guidata handle directly using
sceneG = ieSessionGet('scene window handle')

%% END
