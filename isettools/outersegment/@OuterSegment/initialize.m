% Method to set the default values of all properties
function initialize(obj) 
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    obj.noiseflag = [];
    
    obj.ConeCurrentSignal = [];
    obj.ConeCurrentSignalPlusNoise = [];
    
%     obj.filterKernel = [];
    
%     obj.outputSignal = [];
%     obj.inputSignal = [];
%     obj.time = [];
%     obj.timeConstant = 100;
    
    
end

