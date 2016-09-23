function [uData, hf] = plot(obj, type, varargin)
% Plot function for coneMosaicHex (subclass of coneMosaic)
%
%    [uData, hf] = coneMosaicHex.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure handle or control structure, the meaning of value is%             []            - create plot in new figure
%             figure handle - plot in figure specified by hf
%             'none'        - don't plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
% See coneMosaic.plot for plot types
%
% Example:
%    rgc.mosaic{1}.plot(type)
%
% HJ/BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj');
p.addRequired('type', @isstr);               % Type of plot

p.addParameter('hf', []);                    % figure handle
p.addParameter('oi',[],@isstruct);           % Used for spectral qe

p.parse(obj,type, varargin{:});
hf = p.Results.hf;

uData = [];

%% Figure and axes
if isempty(hf), hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure'), figure(hf);
elseif isgraphics(hf, 'axes'), axes(hf);
end

% set color order so that LMS plots as RGB
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');
    if isgraphics(hf,'axes')
        set(get(hf,'parent'),'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
    else  % Figure
        set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
    end
end

%% Pick up the specialized hex routines here

% If not one of ours, pass on to the base class in the otherwise case.
switch ieParamFormat(type);
    case 'conemosaic'
        % Try the NP plotting routine from here ...
        % It brings the image up in the wrong window because we didn't pass
        % the hf and set the axis.  But the good news is, it gets here and
        % runs.  Next, on to fixing the varargin{} and such.
        plotHexMosaic(obj);  % Default arguments for now
    case 'meanabsorptions'
        disp('NYI')
    case 'movieabsorptions'
        disp('NYI')
    case 'meancurrent'
        disp('NYI')
    case 'moviecurrent'
        disp('NYI')
    otherwise
        % Not one of the hex image types.  So, pass up to the base class
        [uData,hf] = plot@coneMosaic(obj,type,varargin{:});
end
