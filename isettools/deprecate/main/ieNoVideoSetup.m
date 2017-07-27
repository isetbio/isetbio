function ieNoVideoSetup(handles)
% Adjust the vcMainWindow by removing the video button
%
%     ieNoVideoSetup(handles)
%
% As we increase the number of modules, the opening function may have to
% grow. Having the window adjustment take place in a separate routine
% may help us adjust it. 
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('handles')
    error('Requires handles to the vcMainWindow.'); 
end

set(handles.btnVideo,'Visible','off')

u = get(handles.figure1,'units');
set(handles.figure1,'units','normalized');
p = get(handles.figure1,'position');  % .022, .62, .2, .413 are defaults
p(2) = 0.6; p(4) = .35;
set(handles.figure1,'position',p);

p = get(handles.backgroundFrame,'position');
p(2) = p(2) - 0.05;
p(4) = p(4) + 0.03;               %0.213, 0.2, .56, 0.60
set(handles.backgroundFrame,'position',p);

s = 0.08;
p = get(handles.btnScene,'position');
p(2) = p(2)- s; 
set(handles.btnScene,'position',p);

p = get(handles.opticalImage,'position');
p(2) = p(2)- s; 
set(handles.opticalImage,'position',p);

p = get(handles.btnSensorImage,'position');
p(2) = p(2)- s; 
set(handles.btnSensorImage,'position',p);

p = get(handles.btnDisplayImage,'position');
p(2) = p(2)- s; 
set(handles.btnDisplayImage,'position',p);    

set(handles.figure1,'units',u);
end