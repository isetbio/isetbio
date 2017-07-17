function ieMainClose
% Close the all the ISET windows and the ieMainwindow
%
%     ieMainClose
%
% The routine checks for various fields, closes all the main windows
% properly.
%
% Copyright ImagEval Consultants, LLC, 2005.

%% We could probably get rid of vcSESSION in ISETBIO

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

vcSESSION.GUI = [];

%% This closes all the remaining figures in the root.

%  The ISETBIO figures, such as rgcLayer.window and bipolarLayer.window,
%  are closed by this call.
h = findall(0,'Type','Figure');
delete(h);

end
