function obj = osCreate(varargin)
% osCreate: generate an @osLinear, @osBioPhys or @osIdentity object.
% 
% Inputs: none.
% Outputs: the outersegment object, with a value of '0' for the noise flag
% and empty vectors for 'ConeCurrentSignal' and
% 'ConeCurrentSignalPlusNoise'.
% 
% See the initialize method for the @osLinear, @osBioPhys and @osIdentity 
% subclasses for more details of the specific implementations.
%
% 7/2015 JRG


    if nargin == 0
        obj = osLinear();
    elseif nargin == 1
        if strcmpi(varargin{1},'linear');
            obj = osLinear();
        elseif (strcmpi(varargin{1},'biophys') || strcmpi(varargin{1},'rieke'));
            obj = osBioPhys();
        elseif strcmpi(varargin{1},'identity');
            obj = osIdentity();
        else
            obj = osLinear();
        end
    elseif nargin == 3
        if strcmpi(varargin{1},'linear');
            obj = osLinear(varargin{2},varargin{3});
        elseif (strcmpi(varargin{1},'biophys') || strcmpi(varargin{1},'rieke'));
            obj = osBioPhys(varargin{2},varargin{3});
        elseif strcmpi(varargin{1},'identity');
            obj = osIdentity(varargin{2},varargin{3});
        else
            obj = osLinear(varargin{2},varargin{3});
        end
    else        
        obj = osLinear();
    end
    
end

