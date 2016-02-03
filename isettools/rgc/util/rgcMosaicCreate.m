function rgc = rgcMosaicCreate(rgc, varargin)
% rgcCreate: generate an rgcMosaicLinear, rgcMosaicLNP or rgcMosaicGLM object.
%
%  obj = rgcCreate(rgc, mosaicType)
% 


narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% addRequired(p,'rgc',@isa();
addParameter(p,'mosaicType','',@ischar);

p.parse(varargin{:});

%% Spatial RF arrays


switch ieParamFormat(p.Results.mosaicType)
    case{'onparasol'}
        mosaicInd = 1;
    case{'offparasol'}
        mosaicInd = 2;
    case{'onmidget'}
        mosaicInd = 3;
    case{'offmidget'}
        mosaicInd = 4;
    case{'smallbistratified','sbc'}
        mosaicInd = 5;
    otherwise
        mosaicInd = length(rgc.mosaic)+1;
        if mosaicInd == 2 && isempty(rgc.mosaic{1})
            mosaicInd = 1;
        elseif mosaicInd > 5
            mosaicInd = 1;
        end
end
          

% Set key-value pairs
switch ieParamFormat(class(rgc))
    case {'linear','rgclinear'}
        obj = rgcMosaicLinear(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    case {'pool','rgcpool'}
        obj = rgcMosaicPool(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    case {'lnp','rgclnp'}
        obj = rgcMosaicLNP(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    case {'glm','rgcglm'}
        obj = rgcMosaicGLM(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    case {'subunit','rgcsubunit'}
        obj = rgcMosaicSubunit(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    case{'phys','rgcphys'}
        obj = rgcMosaicPhys(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);
    otherwise
        obj = rgcMosaicLinear(rgc, mosaicInd);
        rgcSet(rgc, 'mosaic', obj);        
end

