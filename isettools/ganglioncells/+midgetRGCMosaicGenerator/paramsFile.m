function theParamsFile = paramsFile()
     paramsFilePath = fullfile(isetRootPath, 'ganglioncells');
     paramsFile = 'MidgetRGCMosaicGeneratorParams.mat';
     theParamsFile = fullfile(paramsFilePath, paramsFile);
end