% RGCMosaicConstructor.filepathFor.rawFigurePDFsDir
function theRawFiguresDir = rawFigurePDFsDir()

    p = getpref('isetbio'); 
    if (isempty(p))
        errorString = sprintf('\nYou do not have set your isetbio preferences.\nYou need to do the following.\n(1) Execute ''run(''isetbio/configuration/isetbioLocalHookTemplate'')'' to generate them.\n(2) Edit ''isetbio/isettools/ganglioncells/+RGCMosaicConstructor/+helper/+utils/generateLocalPrefs.m'' and add specific paths for YOUR_COMPUTER_NAME.');
        error(errorString);
    else
        if (~isfield(p, 'rgcResources'))
            errorString = sprintf('\nYou need to add the ''rgcResources'' field in your isetbio prefs.\nEdit ''isetbio/isettools/ganglioncells/+RGCMosaicConstructor/+helper/+utils/generateLocalPrefs.m'' and add specific paths for YOUR_COMPUTER_NAME.');
            error(errorString);
        end
    end

    theRawFiguresDir = p.rgcResources.figurePDFsDir;

 end