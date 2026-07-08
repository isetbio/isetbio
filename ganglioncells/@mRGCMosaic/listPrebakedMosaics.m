function [theRGCmosaicFileNames, prebakedRGCMosaicDir] = listPrebakedMosaics(varargin)
    
    % Parse input
    p = inputParser;
    p.addParameter('type', 'ONcenterLinear', @(x)(ismember(x, {'ONcenterLinear'})));
    p.parse(varargin{:});

 
    switch (p.Results.type)
        case 'ONcenterLinear'
            rgcTypeSubDirectory = 'ONmRGCmosaics';
    end


    prebakedRGCMosaicDir = fullfile('isettools/ganglioncells/data/prebakedRGCmosaics', rgcTypeSubDirectory);
    prebakedRGCMosaicDir = fullfile(isetbioRootPath, prebakedRGCMosaicDir);
    listing = dir(prebakedRGCMosaicDir);

    folders = struct2table(listing);
    theRGCmosaicFileNames = folders(~matches(folders.name,[".","..", ".DS_Store"]),:);

    theRGCmosaicFileNames = theRGCmosaicFileNames{:,1};
end