function figNum = plotSetUpWindow(figNum)
% Initialize position, size and colormap for ISET plot window
%
% Syntax:
%   figNum = plotSetUpWindow([figNum])
%
% Description:
%    This routine could, ultimately, take into account the user's display
%    size or other window positions.
%
% Inputs:
%    figNum - (Optional) The figure number. Default the most recent figure.
%
% Outputs:
%    figNum - The figure that has just been displayed.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: The management of the graphics windows is very poor right now
%      and needs to be corrected. We should allow multiple graph windows
%      and help keep track of them.
%   


% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/11/17  jnm  Formatting
%    01/24/18  jnm  Formatting update to match Wiki

% Examples:
%{
    plotSetUpWindow(1);
%}
%{
    % ETTBSkip - skipping broken examples
    plotSetUpWindow;
%}


if notDefined('figNum')
    figNum = vcSelectFigure('GRAPHWIN'); 
    set(figNum, 'Units', 'Normalized', 'Position', ...
        [0.5769, 0.0308, 0.4200, 0.4200]);
    set(figNum, 'Name', 'ISET GraphWin', 'NumberTitle', 'off');
else
    figure(figNum); 
end

colormap('default')
clf

return;