function figHdl = vcNewGraphWin(figHdl, fType, varargin)
% Open a window for plotting
%
% Syntax:
%   figHdl = vcNewGraphWin([figHdl], [fType], [varargin])
%
% Description:
%    A graph window figure handle is returned and stored in the currernt
%    vcSESSION.GRAPHWIN entry.
%
%    The varargin is a set of (param, val) pairs that can be used for
%    set(gcf, param, val);
%
%    A few figure shapes can be defaulted
%      fType:  Default - Matlab normal figure position
%              upper left    Simple
%              tall          (for 2x1 format)
%              wide          (for 1x2 format)
%              upperleftbig  (for 2x2 format)
%      This list may grow.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcNewGraphWin.m' into the Command Window.
%
% Inputs:
%    figHdl   - (Optional) Handle. The figure handle. Default calls figure.
%    fType    - (Optional) String. The figure shape. Supported options
%               include the following: 
%        'Default'       Matlab normal figure position
%        'upper left'    Simple. (Actual) Default.
%        'tall'          (for 2x1 format)
%        'wide'          (for 1x2 format)
%        'upperleftbig'  (for 2x2 format)
%    varargin - (Optional) Additional variable(s) as necessary to perform
%               the function call. These can be seen in the key/values
%               section below.
%
% Outputs:
%    figHdl   - Handle. The figure handle.
%
% Optional key/value pairs:
%    **NEEDS TO BE FILLED OUT**
%
% Notes:
%    * TODO: Add more default position options
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/11/18  jnm  Formatting

% Examples:
%{
    vcNewGraphWin;
    vcNewGraphWin([], 'upper left')
    vcNewGraphWin([], 'tall')
    vcNewGraphWin([], 'wide')
    vcNewGraphWin([], 'upper left big')

    vcNewGraphWin([], 'wide', 'Color', [0.5 0.5 0.5])
%}

if notDefined('figHdl'), figHdl = figure; end
if notDefined('fType'), fType = 'upper left'; end

set(figHdl, 'Name', 'ISET GraphWin', 'NumberTitle', 'off');
% set(figHdl, 'CloseRequestFcn', 'ieCloseRequestFcn');
set(figHdl, 'Color', [1 1 1]);

% Position the figure
fType = ieParamFormat(fType);
switch(fType)
    case 'upperleft'
        set(figHdl, 'Units', 'normalized', ...
            'Position', [0.007 0.55  0.28 0.36]);
    case 'upperright'
        set(figHdl, 'Units', 'normalized', ...
            'Position', [0.30 0.55   0.28 0.36]);
    case 'tall'
        set(figHdl, 'Units', 'normalized', ...
            'Position', [0.007 0.055 0.28 0.85]);
    case 'wide'
        set(figHdl, 'Units', 'normalized', ...
            'Position', [0.007 0.62  0.7  0.3]);
    case 'upperleftbig'
        % Like upperleft but bigger
        set(figHdl, 'Units', 'normalized', ...
            'Position', [0.007 0.40  0.40 0.50]);

    otherwise % Matlab default
end

%% Parse the varargin arguments
if ~isempty(varargin)
    n = length(varargin);
    if ~mod(n, 2)
        for ii = 1:2:(n - 1)
            set(figHdl, varargin{ii}, varargin{ii + 1});
        end
    end
end

%% Store some information. Not sure it is needed; not much used.
% ieSessionSet('graphwinfigure', figHdl);
% ieSessionSet('graphwinhandle', guidata(figHdl));

end