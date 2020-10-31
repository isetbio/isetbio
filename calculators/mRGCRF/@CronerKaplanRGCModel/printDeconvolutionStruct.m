function printDeconvolutionStruct(deconvolutionStruct)
    
    deconvolutionStruct{1,1}.metaData
    coneInputConfigNames = keys(deconvolutionStruct{1,1}.data);
    
    for coneInputConfigIndex = 1:numel(coneInputConfigNames)
        sprintf('config: %s\n', coneInputConfigNames{coneInputConfigIndex});
        deconvolutionStruct{1,1}.data(coneInputConfigNames{coneInputConfigIndex})
    end
    
end

