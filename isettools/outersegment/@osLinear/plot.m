function h = plot(os, pType, varargin)
% Plots the osLinear properties and then sends on plot@outersegment
%
% Inputs: 
%   os - osLinear object
% 
% pType for osLinear only:
%  current filters - Converts absorptions to photo current for linear model
%
% Outputs: 
%    h is a handle to the plot window
%
% Properties that can be plotted:
%
% Examples:
%   osL.plot(absorptions,'type','isomerizations')
%   osL.plot(absorptions,'type','current')
%   osL.plot(absorptions,'type','filter kernels')
%   osL.plot(absorptions,'type','all')
%
% BW (c) isetbio team, 2016

%% Check for the number of arguments and create parser object.

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('os',@(x)(isa(os,'osLinear')));
p.addRequired('pType',@ischar);

p.addParameter('meancurrent',[],@isvector);

p.parse(os, pType, varargin{:});
meancurrent = p.Results.meancurrent;  % Background current for linear model

% No additional parameter/values yet
% Could put in for one cone class

%% Choosing the plot type or send on to base class

switch ieParamFormat(pType)
        
    case {'currentfilters'}
        % Plot linear temporal filters for L, M and S cones.

        h = vcNewGraphWin;
        tSamples = (1:size(os.lmsConeFilter,1))*os.timeStep*1000;
        
        plot(tSamples,os.lmsConeFilter(:,1),'r-', ...
            tSamples,os.lmsConeFilter(:,2),'g-', ...
            tSamples,os.lmsConeFilter(:,3),'b-');
        xlabel('Time (ms)'); ylabel('Current (pA)');
        grid on;
        l = cell(1,3);
        if ~isempty(meancurrent)
            c = {'L','M','S'};
            for ii=1:3
                l{ii} = sprintf('%s cone (mean %.2f pA)',c{ii},meancurrent(ii));
            end
        else
            l = {'L cone','M cone','S cone'};
        end
        legend(l);
        title('Absorption impulse response');
        
    otherwise
        plot@outerSegment(os,pType);
end

end





