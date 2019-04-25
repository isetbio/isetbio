function [uData, h] = plot(os, pType, varargin)
% OSLINEAR.PLOT - Plot osLinear properties or sends to @outersegment.plot
%
% Syntax:
%   [uData, h] = plot(os, pType, [varargin])
%
% Description:
%    Plot the osLinear properties or send them to @outersegment.plot.
%
%    Examples are contained in the code. To access, type 'edit plot.m'
%    into the Command Window.
%
% Inputs:
%	 os       - osLinear object
%    pType    - Plot Type. 'current filters' is a type specific to
%               osLinear, which converts absorptions to photo current for a
%               linear model. Other options include 'isomerizations',
%               'current', 'filter kernels', 'all'.
%    varargin - (Optional) Additional information for the plot.
%
% Outputs:
%    uData    - User data for values in the plot
%               (use uData = get(gca, 'userdata);)
%    h        - Handle to the plot window
%
% Optional key/value pairs:
%    'cmosaic'      - Parent @coneMosaic object
%    'meancurrent'  - Background current for linear model
%

% History:
%    xx/xx/16  BW   (c) isetbio team, 2016
%    02/13/18  jnm  Formatting
%    04/07/18  dhb  Skip broken example.

% Examples:
%{
    % ETTBSkip. Even when you initialize the osL object, the functions do
    % not work.
    osL = osCreate;
    osL.plot(absorptions, 'type', 'isomerizations')
    osL.plot(absorptions, 'type', 'current')
    osL.plot(absorptions, 'type', 'filter kernels')
    osL.plot(absorptions, 'type', 'all')
%}

%% Check for the number of arguments and create parser object.

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('os', @(x)(isa(os, 'osLinear')));
p.addRequired('pType', @ischar);

% Sometime we need these
p.addParameter('cmosaic', [], @(x)(isa(x, 'coneMosaic')));
p.addParameter('meancurrent', [], @isvector);

p.parse(os, pType, varargin{:});
meancurrent = p.Results.meancurrent;  % Background current for linear model
cmosaic = p.Results.cmosaic;

% No additional parameter/values yet
% Could put in for one cone class

%% Choosing the plot type or send on to base class
switch ieParamFormat(pType)
    case {'impulseresponse', 'currentfilters'}
        % Plot linear temporal filters for L, M and S cones.
        h = vcNewGraphWin;

        if isempty(os.lmsConeFilter), os.linearFilters(cmosaic); end
        tSamples = os.timeAxis;
        uData.t = tSamples;
        uData.y = os.lmsConeFilter;

        plot(tSamples, os.lmsConeFilter(:, 1), 'r-', ...
            tSamples, os.lmsConeFilter(:, 2), 'g-', ...
            tSamples, os.lmsConeFilter(:, 3), 'b-');
        xlabel('Time (sec)');
        ylabel('Current (pA)');
        grid on;

        % Create legend and title
        l = cell(1, 3);
        if ~isempty(meancurrent)
            c = {'L', 'M', 'S'};
            for ii = 1:3
                l{ii} = sprintf('%s cone (mean %.2f pA)', c{ii}, ...
                    meancurrent(ii));
            end
        else
            l = {'L cone', 'M cone', 'S cone'};
        end
        legend(l);
        title('Absorption impulse response');

    otherwise
        [uData, h] = plot@outerSegment(os, pType);
end

set(gca, 'userdata', uData);

end
