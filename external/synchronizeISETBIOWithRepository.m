% Function to synchronize external functions in isetbio with their origins. 
% 
% Examples are provided in the code.
%
%  11/20/2014 npc   Wrote it.
%  4/27/2015 npc    Modification to support synchronization of any repository with isetbio, not just PTB 
%

% Examples:
%{
    % ETTBSkip
    % Running these can get you stale files if your source
    % is out of date, so don't auto-exacute these examples.

    % Synchronize ISETBIO's external PTB routines
    synchronizeISETBIOWithRepository('PTB_DHB');

    % Synchronize ISETBIO's external RTB-3 routines
    % THIS SHOULD BE UPDATED TO RTB-4
    synchronizeISETBIOWithRepository('RT3_DHB');

    % Synchronize ISETBIO's external BrainardLab Toolbox routines
    synchronizeISETBIOWithRepository('BLTB_DHB');
%}
function synchronizeISETBIOWithRepository(repositoryName)
   
    if (nargin ~= 1)
        error('Usage: synchronizeISETBIOWithRepository(''repositoryName''), i.e., synchronizeISETBIOWithRepository(''PTB'').');
    end

    switch repositoryName
        case 'PTB_DHB'
            srcDir = '/Users/dhb/Documents/MATLAB/toolboxes/Psychtoolbox-3/Psychtoolbox';
            dstDir = '/Users/dhb/Documents/MATLAB/toolboxes/isetbio/external/psychtoolbox';
        case 'RT3_DHB'
            srcDir = '/Users/Shared/Matlab/Toolboxes/RenderToolbox3';
            dstDir = '/Users/dhb/Documents/MATLAB/toolboxes/isetbio/external/rendertoolbox3';
        case 'BLTB_DHB'
            srcDir = '/Users/dhb/Documents/MATLAB/toolboxes/BrainardLabToolbox';
            dstDir = '/Users/dhb/Documents/MATLAB/toolboxes/isetbio/external/brainardlabtoolbox';
        otherwise
            error('No information about how to synchronize ''%s''. Please update ''synchronizeISETBIOWithRepository.m''.', repositoryName);
    end

    % Check that directories exist
    if (~exist(srcDir, 'dir'))
        error('Source directory, ''%s'', does not exist.', srcDir);
    end

    if (~exist(dstDir, 'dir'))
        error('Destination directory, ''%s'', does not exist.', dstDir);
    end
    
    % Get names of sub-directories of $isetbio/external/ptb
    fprintf('\nChecking contents of destination directory (%s) to determine which files to sync.\n', dstDir);
    dirContents = dir(dstDir);
    ignoredDirs = {'+ptb', '', '.', '..'};
    destDirs = {}; srcDirs = {};
    for k = 1:numel(dirContents)
        if (dirContents(k).isdir)
            if ismember(dirContents(k).name, ignoredDirs)
                continue;
            end
            fullDirName = sprintf('%s/%s', dstDir, dirContents(k).name);
            destDirs{numel(destDirs)+1} = fullDirName;
            fullDirName = sprintf('%s/%s', srcDir, dirContents(k).name);
            srcDirs{numel(srcDirs)+1} = fullDirName;
        else
            fprintf('%s is not a directory. Will not be synchronized.\n', dirContents(k).name);
        end
    end
    
    % Recursively scan all directories of $isetbio/external/ptb and get all filenames
    filesToSynchronize.src = {}; filesToSynchronize.dest = {};
    for k = 1:numel(destDirs)
        filesToSynchronize = updateFilesToSynchronize(destDirs{k}, srcDirs{k}, filesToSynchronize);
    end
    
    % Rest counters
    MFiles.different= 0;     MFiles.same     = 0;
    MATFiles.different = 0;  MATFiles.same   = 0;
    TXTFiles.different = 0;  TXTFiles.same = 0;
    DOCXFiles.different = 0; DOCXFiles.same = 0;
    UNKNOWNFORMATFiles.different = 0; UNKNOWNFORMATFiles.same = 0;

    % Do the synchronization
    fprintf('\nPlease wait. Diffing %d files ...\n', numel(filesToSynchronize.dest));
    report = {};
    for k = 1:numel(filesToSynchronize.dest)
        destinationFile = filesToSynchronize.dest{k};
        sourceFile      = filesToSynchronize.src{k};
        if (exist(sourceFile, 'file') ~= 2) && (isempty(strfind(sourceFile, '.DS_Store')))
            error('%s does not exist !!', sourceFile);
        else 
            [content_differs, printout] = system(sprintf('diff --side-by-side  %s %s', sourceFile, destinationFile));
            [pathstr,name,ext] = fileparts(sourceFile);
            if (strcmp(ext, '.m'))
                %fprintf('Found source m file %s.\n', sourceFile);
                if (content_differs ~= 0)
                    MFiles.different = MFiles.different + 1;
                else
                    MFiles.same = MFiles.same + 1;
                end
            elseif (strcmp(ext, '.mat'))
                %fprintf('Found source mat file %s.\n', sourceFile);
                if (content_differs ~= 0)
                    MATFiles.different = MATFiles.different + 1;
                else
                    MATFiles.same = MATFiles.same + 1;
                end
            elseif (strcmp(ext, '.txt'))
                %fprintf('Found source txt file %s.\n', sourceFile);
                if (content_differs ~= 0)
                    TXTFiles.different = TXTFiles.different + 1;
                else
                    TXTFiles.same = TXTFiles.same + 1;
                end
            elseif (strcmp(ext, '.docx'))
                %fprintf('Found source docx file %s.\n', sourceFile);
                if (content_differs ~= 0)
                    DOCXFiles.different = DOCXFiles.different + 1;
                else
                    DOCXFiles.same = DOCXFiles.same + 1;
                end
            else
                %fprintf(2,'Found source file (unknown format) %s.\n', sourceFile);
                if (content_differs ~= 0)
                    UNKNOWNFORMATFiles.different = UNKNOWNFORMATFiles.different + 1;
                else
                   UNKNOWNFORMATFiles.same =  UNKNOWNFORMATFiles.same + 1;
                end
            end
            
            if (content_differs ~= 0)
                system(sprintf('cp %s %s', sourceFile, destinationFile));
                report{numel(report)+1} = sprintf('%4d. %-150s %s %s\n', numel(report)+1, destinationFile, '<-', sourceFile);
            end
        end
    end
    fprintf(' Done !\n');
     
    fprintf('\nScaned %4d total files.', numel(filesToSynchronize.dest));
    fprintf('\nOut of %4d M-files            %d needed to be updated.', MFiles.different    + MFiles.same,    MFiles.different);
    fprintf('\nOut of %4d MAT-files          %d needed to be updated.', MATFiles.different  + MATFiles.same,  MATFiles.different);
    fprintf('\nOut of %4d TXT-files          %d needed to be updated.', TXTFiles.different  + TXTFiles.same,  TXTFiles.different);
    fprintf('\nOut of %4d DOCX-files         %d needed to be updated.', DOCXFiles.different + DOCXFiles.same, DOCXFiles.different);
    fprintf('\nOut of %4d other-format-files %d needed to be updated.', UNKNOWNFORMATFiles.different + UNKNOWNFORMATFiles.same, UNKNOWNFORMATFiles.different);        
    fprintf('\n');
    
    fprintf('\nUpdated files:\n');
    for k = 1:numel(report)
        fprintf('%s', char(report{k}));
    end
end

    
% Method to recursively scan directory destDirName and get all file names
function filesToSynchronize = updateFilesToSynchronize(destDirName, srcDirName, filesToSynchronize)  
    % get contents of dir name
    dirContents = dir(destDirName);
    ignoredDirs = {'.', '..'};
    
    for k = 1:numel(dirContents)
        if (dirContents(k).isdir)
            if ismember(dirContents(k).name, ignoredDirs)
                continue;
            end
            fullDestDirName = sprintf('%s/%s', destDirName, dirContents(k).name);
            fullSrcDirName  = sprintf('%s/%s', srcDirName, dirContents(k).name);
            filesToSynchronize = updateFilesToSynchronize(fullDestDirName, fullSrcDirName, filesToSynchronize);
        else
            fullDestFileName = sprintf('%s/%s',destDirName, dirContents(k).name);
            fullSrcFileName  = sprintf('%s/%s',srcDirName, dirContents(k).name);
            filesToSynchronize.dest{numel(filesToSynchronize.dest)+1} = fullDestFileName;
            filesToSynchronize.src{numel(filesToSynchronize.src)+1} = fullSrcFileName;
        end
    end
end