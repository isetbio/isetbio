function ieCloseRequestFcn
% Deprecated:  Method used for closing ISET graph windows
%
% Syntax:
%   ieCloseRequestFcn
%
% Description:
%    When we close a graph window, we try to clean up the vcSESSION
%    information whenever we close the currently selected graph win.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/28/18  jnm  Formatting

curFig = gcf;
% curGraphWin = ieSessionGet('graphwinfigure');

% If this is the current graph window, then we set the figure information
% to empty before exiting.
% if (curFig == curGraphWin)
%     ieSessionSet('graphwinfigure',[]);
%     ieSessionSet('graphwinhandle',[]);
% end

closereq;

return;