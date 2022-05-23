% Controlling window properties in Matlab.
%
% Description:
%    It is possible to change some of the GUI appearance properties from
%    the command line.  Here we illustrate how to make windows invisible, 
%    change the white point, and change the font size.
%
% Notes:
%    1) This might be elaborated a bit, to show explicitly more things that
%       you might want to do.
%

% History:
%    XX/XX/10       Copyright ImagEval Consultants, LLC, 2010
%    11/22/18  JNM  Formatting

%% Initialize
ieInit;

%% Opening and closing windows from the command line
% Create an isetbio window
scene = sceneCreate;
sceneW = sceneWindow(scene);

sceneW.Visible = 'off';
pause(1);
sceneW.Visible = 'on';
pause(1);


%% Interacting with the scene handles
% It is also possible to make adjustments to the display by interacting
% with the Matlab handle graphics.  To get the handle to the scene figure,
% you can run
sceneF = ieSessionGet('scene figure');

% The variable sceneF is the handle to the figure
get(sceneF)

% The guidata are available here
guidata(sceneF)

% Or you can get the guidata handle directly using
sceneG = ieSessionGet('scene window handle');

%% END
