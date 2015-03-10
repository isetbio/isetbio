function ieCloseRequestFcn
%Method used for closing ISET windows
%
%   ieCloseRequestFcn
%
% Purpose:
%   When we close a graph window, we try to clean up the vcSESSION
%   information whenever we close the currently selected graph win.
%
% Copyright ImagEval Consultants, LLC, 2003.

curFig = gcf;
curGraphWin = ieSessionGet('graphwinfigure');

% If this is the current graph window, then we set the figure information
% to empty before exiting.
if (curFig == curGraphWin)
    ieSessionSet('graphwinfigure',[]);
    ieSessionSet('graphwinhandle',[]);
end

closereq;

return;