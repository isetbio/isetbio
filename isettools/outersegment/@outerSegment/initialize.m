% Method to set the default values of all properties
function initialize(obj) 
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    obj.noiseflag = [];
    
    obj.ConeCurrentSignal = [];
    obj.ConeCurrentSignalPlusNoise = [];
    
%     filterKernel(obj);
    
end

