function initialize(obj)
% initialize: an abstract method of @outerSegment. The initialize methods 
% of @osLinear and @osBioPhys set the default values of object of those
% subclasses. 
% 
% Inputs: the outersegment object.
% Outputs: the outersegment object, with a value of '0' for the noise flag
% and empty vectors for 'ConeCurrentSignal' and
% 'ConeCurrentSignalPlusNoise'.
% 
% See the initialize method for the @osLinear and @osBioPhys subclasses for
% more details of the specific implementations.
%
% 7/2015 JRG

    obj.noiseFlag = 0;
    
    obj.ConeCurrentSignal = [];
    obj.ConeCurrentSignalPlusNoise = [];
    
end

