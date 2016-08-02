function obj = mosaicSet(obj, param, val, varargin)
%  Sets a property for an rgcLNP object.
% 
%   rgc.mosaic = mosaicSet(rgc.mosaic, param, value, varargin)
%  
% Inputs: 
% 
%   obj    - rgc object
%   param  - parameter string
%   val    - parameter value
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%    obj with property set appropriately
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated

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
        % The baseline rate of the likelihood of spiking, which creates a
        % nonzero firing rate in response to no input
        for ci1 = 1:size(obj.tonicDrive,1)
            for ci2 = 1:size(obj.tonicDrive,2)
                obj.tonicDrive{ci1,ci2} = params.value;
            end
        end
    case{'generatorfunction'}
        % Converts the conditional intensity to the likelihood of observing
        % a spike within a given time bin, used also in GLM case.
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

