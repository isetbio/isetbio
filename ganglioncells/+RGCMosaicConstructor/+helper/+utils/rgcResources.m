%
% RGCMosaicConstructor.helper.utils.rgcResources
%
function theRGCresources = rgcResources()
    p = getpref('isetbio'); 
    if (isempty(p))
        errorString = sprintf('\nYou do not have set your isetbio preferences.\nYou need to do the following.\n(1) Execute ''run(''isetbio/configuration/isetbioLocalHookTemplate'')'' to generate them.\n(2) Edit ''isetbio/ganglioncells/+RGCMosaicConstructor/+helper/+utils/generateLocalPrefs.m'' and add specific paths for YOUR_COMPUTER_NAME.');
        error(errorString);
    end
        
    if (~isfield(p, 'rgcResources'))
        errorString = sprintf('\nYou need to add the ''rgcResources'' field in your isetbio prefs.\nEdit ''isetbio/ganglioncells/+RGCMosaicConstructor/+helper/+utils/generateLocalPrefs.m'' and add specific paths for YOUR_COMPUTER_NAME.');
        error(errorString);
    end

    theRGCresources = p.rgcResources;

    % Output directories must be absolute. A relative value (such as the
    % 'FullPathTo...' placeholder that generateLocalPrefs used to default to)
    % resolves against whatever the current directory happens to be, which for
    % a tutorial launched with run() is the source tree itself. Fall back to
    % isetbio's gitignored local/ directory rather than writing there.
    theRGCresources.intermediateDataDir = localResolveOutputDir( ...
        theRGCresources, 'intermediateDataDir', ...
        fullfile(isetbioRootPath, 'local', 'rgcResources', 'intermediateFiles'));

    theRGCresources.figurePDFsDir = localResolveOutputDir( ...
        theRGCresources, 'figurePDFsDir', ...
        fullfile(isetbioRootPath, 'local', 'rgcResources', 'figures'));
end

function theDir = localResolveOutputDir(theRGCresources, fieldName, defaultDir)
    theDir = '';
    if (isfield(theRGCresources, fieldName))
        theDir = theRGCresources.(fieldName);
    end

    if (localIsAbsolutePath(theDir)), return; end

    persistent alreadyWarned
    if (isempty(alreadyWarned)), alreadyWarned = struct(); end
    if (~isfield(alreadyWarned, fieldName))
        alreadyWarned.(fieldName) = true;
        fprintf(['\nisetbio prefs give a non-absolute ''%s'' (''%s'').\n' ...
                 'Writing to ''%s'' instead. Run ' ...
                 'RGCMosaicConstructor.helper.utils.generateLocalPrefs() to set your own.\n\n'], ...
                 fieldName, theDir, defaultDir);
    end

    theDir = defaultDir;
end

function tf = localIsAbsolutePath(thePath)
    if (isempty(thePath)) || (~(ischar(thePath) || isStringScalar(thePath)))
        tf = false;
        return;
    end
    if (ispc)
        tf = ~isempty(regexp(thePath, '^([A-Za-z]:[\\/]|\\\\)', 'once'));
    else
        tf = startsWith(thePath, '/') || startsWith(thePath, '~');
    end
end
