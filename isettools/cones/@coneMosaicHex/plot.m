function [uData, hf] = plot(obj, type, varargin)
% Plot function for coneMosaicHex (subclass of coneMosaic)
%
% Syntax:
%	[uData, hf] = coneMosaicHex.plot(type, [varargin])
%   [uData, hf] = plot(obj, type, [varargin])
%
% Inputs:
%    obj  - The cone mosaic hex object
%    type - A string denoting the type of the plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
% Optional key/value pairs:
%   'hf'  - figure handle or control structure, the meaning of value is
%               []            - create plot in new figure
%               figure handle - plot in figure specified by hf
%               'none'        - don't plot
%
% See Also:
%    coneMosaic.plot for plot types
%

% History:
%    xx/xx/16  HJ/BW  ISETBIO TEAM, 2016
%    02/20/18  jnm    Formatting
%    04/07/18  dhb    Skip broken example.

% Examples:
%{
    % ETTBSkip. This needs a mosaic defined before it could possibly work.
    rgc.mosaic{1}.plot(type)
%}

%% parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj');
p.addRequired('type', @isstr);        % Type of plot

p.addParameter('hf', []);             % figure handle
p.addParameter('oi', [], @isstruct);  % Used for spectral qe

p.parse(obj, type, varargin{:});
hf = p.Results.hf;

uData = [];

%% Figure and axes
if isempty(hf)
    hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure')
    figure(hf);
elseif isgraphics(hf, 'axes')
    axes(hf);
end

% Set color order so that LMS plots as RGB
% Matlab default is 7 colors, and we reorder
% If the user has changed the default, we leave it alone.
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');  
    if size(co,1) == 7  && ~isequal(co(1,:),[1,0,0])  % Figure        
        if isgraphics(hf, 'axes')
            set(get(hf, 'parent'), ...
                'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
        else
            set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
        end
    end
end

%% Pick up the specialized hex routines here
% If not one of ours, pass on to the base class in the otherwise case.
switch ieParamFormat(type)
    case 'conemosaic'
        % Try the NP plotting routine from here ...
        % It brings the image up in the wrong window because we didn't pass
        % the hf and set the axis. But the good news is, it gets here and
        % runs. Next, on to fixing the varargin{} and such.
        %
        % [Note: DHB - I'm not sure I had to change this, but passing
        % varargin{:} onwards is a bit confusing to me, as hf can get
        % passed to the next routine twice if you're not careful. Probably
        % should explicitly set anything that we want plotHexMosaic to
        % inherit here.]

        % Default arguments for now
        % plotHexMosaic(obj, 'hf', hf, varargin{:});
        plotHexMosaic(obj, 'hf', hf);

    case 'meanabsorptions'
        if isempty(obj.absorptions)
            plotHexMosaic(obj,'hf',hf);
        else
            plotHexMeanImage(obj, 'type', 'absorptions');
        end

    case 'movieabsorptions'
        uData = movieHex(obj, 'type', 'absorptions');


    case 'meancurrent'
        plotHexMeanImage(obj, 'type', 'current');

    case 'moviecurrent'
        uData = movieHex(obj, 'type', 'current');


    otherwise
        % Not one of the hex image types. So, pass up to the base class
        [uData, hf] = plot@coneMosaic(obj, type, varargin{:});
end
