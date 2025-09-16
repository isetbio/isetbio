function theAugmentedFilePath = augmentedPathWithSubdirs(theRootDir, filePath, varargin)

    p = inputParser;
    p.addParameter('generateMissingSubDirs', false, @islogical);
    % Execute the parser
    p.parse(varargin{:});
    generateMissingSubDirs = p.Results.generateMissingSubDirs;

    [augmentedFullPath, theFileName, ext] = fileparts(filePath);
    theFileName = sprintf('%s%s',theFileName, ext);
    augmentedSubDirs = split(augmentedFullPath, filesep);

    if (isempty(augmentedSubDirs))
        theAugmentedFilePath = fullfile(theRootDir, theFileName);
    else
        augmentedPath = theRootDir;
        for iSubDir = 1:numel(augmentedSubDirs)
            if (numel(augmentedSubDirs{iSubDir})==0)
                continue;
            end
            augmentedPath = fullfile(augmentedPath, augmentedSubDirs{iSubDir});
            if (~isfolder(augmentedPath))
                if (~generateMissingSubDirs)
                    error('File ''%s'', does not exist. ', augmentedPath);
                end
                fprintf(2,'Will generate subdir ''%s'' within %s\n', augmentedSubDirs{iSubDir}, strrep(augmentedPath, augmentedSubDirs{iSubDir}, ''));
                fprintf(2,'Hit enter to continue...');
                pause;
                mkdir(augmentedPath);
            end
        end % for iSubDir
        theAugmentedFilePath = fullfile(augmentedPath, theFileName);
    end

end
