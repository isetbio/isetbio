function obj = mosaicSet(obj, param, val, varargin)
% mosaicSet for LNP subclass, superclass is @rgcMosaic
% 
%   rgc.mosaic = mosaicSet(rgc.mosaic, param, value, varargin)
%  
% Inputs: 
%   rgc object, key-value pair of property and value to which it is
%   being set.
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% ISETBIO Team, 2016

%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.addRequired('param',@ischar);
p.addRequired('val');

% Parse and put results into structure p.
p.parse(param,val, varargin{:});
param = ieParamFormat(p.Results.param);
val   = p.Results.val;

%% Set key-value pairs
switch param
    
    % LNP subclass parameters
    case{'tonicdrive'}
        %
        for ci1 = 1:size(obj.tonicDrive,1)
            for ci2 = 1:size(obj.tonicDrive,2)
                obj.tonicDrive{ci1,ci2} = params.value;
            end
        end
    case{'generatorfunction'}
        % @JRG:  
        % Does this convert linear to nonlinear for GLM and also LNP?
        obj.generatorFunction = params.value;   
    case{'postspikefilter'}
        obj.postSpikeFilter = params.value;
    case{'responsevoltage'}
        % The nonlinear voltage response after application of the generator
        % function and the spike coupling responses is represented here
        obj.responseVoltage = params.value;
        
    otherwise
        % Super class parameters
        mosaicSet@rgcMosaic(obj,param,val,varargin{:});
        
end

