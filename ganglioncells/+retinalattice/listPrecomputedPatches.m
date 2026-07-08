function [theLatticePatchFileNames, latticeGalleryDir] = listPrecomputedPatches(varargin)

    % Parse input
    p = inputParser;
    p.addParameter('withGenerationProgressHistory', false, @islogical);
    p.parse(varargin{:});

    withGenerationProgressHistory = p.Results.withGenerationProgressHistory;

    latticeGalleryDir = fullfile(isetbioRootPath, 'isettools/ganglioncells/data/lattices');
    listing = dir(latticeGalleryDir);

    folders = struct2table(listing);
    allLatticePatchFileNames = folders(~matches(folders.name,[".","..", ".DS_Store"]),:);

    allLatticePatchFileNames = allLatticePatchFileNames{:,1};
    theLatticePatchFileNames = cell(0,0);

    if (withGenerationProgressHistory)
        for iFile = 1:numel(allLatticePatchFileNames)
            theFileName = allLatticePatchFileNames{iFile};
            if (contains(theFileName, 'progress'))
                theLatticePatchFileNames{size(theLatticePatchFileNames,2)+1} = theFileName;
            end
        end
    else
       for iFile = 1:numel(allLatticePatchFileNames)
            theFileName = allLatticePatchFileNames{iFile};    
            if (~contains(theFileName, 'progress'))
                theLatticePatchFileNames{size(theLatticePatchFileNames,2)+1} = theFileName;
            end
        end
    end

    theLatticePatchFileNames = theLatticePatchFileNames';
end
