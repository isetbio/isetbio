function val = osGet(obj, param)
% Gets isetbio outersegment object parameters.
% 
% Parameters:
%       {'noiseFlag'} -  sets current as noise-free ('0') or noisy ('1')
%       {'sConeFilter'} - the linear filter for S-cone temporal response
%       {'mConeFilter'} - the linear filter for M-cone temporal response
%       {'lConeFilter'} - the linear filter for L-cone temporal response
%       {'patchSize'} - cone current as a function of time
%       {'timeStep'} - noisy cone current signal
%       {'size'} - array size of photon rate
%       {'coneCurrentSignal'} - cone current as a function of time
% 
% osGet(adaptedOS, 'noiseFlag')
% 
% 8/2015 JRG NC DHB


% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 
% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
% narginchk(0, Inf);
% p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;
% 
% % Make key properties that can be set required arguments, and require
% % values along with key names.
% allowableFieldsToSet = {...
%     'noiseflag',...
%     'sconefilter',...
%     'mconefilter',...
%     'lconefilter',...
%     'conespacing',...
%     'conesampling',...
%     'conecurrentsignal'};
% p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
% 
% % Define what units are allowable.
% % allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% % p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));
% % p.addParameter('sconefilter',0,@isnumeric);
% % p.addParameter('mconefilter',0,@isnumeric);
% % p.addParameter('lconefilter',0,@isnumeric);
% 
% % Parse and put results into structure p.
% p.parse(varargin{:}); params = p.Results;

switch ieParamFormat(param)

    case{'sconefilter'}
        val = obj.sConeFilter;
        
    case{'mconefilter'}
        val = obj.mConeFilter;
    
    case{'lconefilter'}
        val = obj.lConeFilter;        
        
    case{'patchsize'}
        val = obj.patchSize;
        
    case{'timestep'}
        val = obj.timeStep;
        
    case{'size'}
        val = size(obj.coneCurrentSignal);
        
    case{'conecurrentsignal','current'}
        val = obj.coneCurrentSignal;
end

