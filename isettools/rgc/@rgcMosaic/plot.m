function [uData, hf] = plot(obj, type, varargin)
% Plot function for rgcMosaic
%    [uData, hf] = rgcMosaic.plot(type, varargin)
%
% Required input
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure or axis handle 
%      By default, this is the figure handle of the rgcMosaic window
%      It can also be one of the two axes in the window.
%      Otherwise
%             []            - create plot in new figure
%             figure handle - an alternative figure
%             'none'        - don't plot, only generate uData
%
% Outputs:
%   uData - Computed user data
%   hf    - figure handle
%
% Plot type can be chosen from
%   'rgc mosaic'          - Color image of the cone arrangement
%   'spike mean image'     - Gray scale movie of current
%   'spike movie'   - Cone photocurrent graphs
%
% Example:
%
% BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

% What about obj?
p.addRequired('type', @isstr);                        % Type of plot
p.addParameter('hf', obj.figureHandle, @isgraphics);  % figure handle

p.parse(type, varargin{:});
hf = p.Results.hf;

uData = [];

% Select place for plot to appear
% The mosaic window has a response image and a geometry image
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

fprintf('Plot for %s class\n',class(obj));

switch ieParamFormat(type)
    case 'spikemeanimage'
        disp('spike image')
    case 'spikemovie'
        disp('Spike Movie')
    case 'psth'
        disp('PSTH')
    otherwise
end

end
