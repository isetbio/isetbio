function ieMainClose
% Close the all the ISET windows and the ieMainwindow
%
%     ieMainClose
%
% The routine checks for various fields, closes all the main windows
% properly.
%
% Copyright ImagEval Consultants, LLC, 2005.

global vcSESSION

if ~checkfields(vcSESSION,'GUI'); closereq; return; end

if checkfields(vcSESSION.GUI,'vcSceneWindow','hObject')
    sceneWindow;
    sceneClose;
end

if checkfields(vcSESSION.GUI,'vcOptImgWindow','hObject')
    oiWindow;
    oiClose;
end

% We need to try this, or perhaps to store the handles in vcSESSION.
% findall(), findobj()

vcSESSION.GUI = [];

% This closes the bipolarLayer and rgcLayer windows, for some reason
% If we could find out how many there are, we could call this again and
% again.  But the possible code above doesn't return any of the windows!
% Why? (BW)
closereq;   

return;