function  [pdfExportSubDir, theSummaryPDFfileName] = optimizationResultsSummaryPDF(theOptimizationResultsFileName, varargin)

    p = inputParser;
    p.addParameter('extraSubDirPath', '', @ischar);
    p.addParameter('generateMissingSubDirs', false, @islogical);
    % Execute the parser
    p.parse(varargin{:});
    extraSubDirPath = p.Results.extraSubDirPath;
    generateMissingSubDirs = p.Results.generateMissingSubDirs;

    [pdfExportSubDir, theSummaryPDFfileName] = fileparts(theOptimizationResultsFileName);

    theOptimizationResultsFileName = strrep(theOptimizationResultsFileName, pdfExportSubDir, '');
    theSummaryPDFfileName = sprintf('%s.pdf', theSummaryPDFfileName);
    theOptimizationResultsFileName = fullfile('optResultsSummaryPDFs',theOptimizationResultsFileName);
    
    if (contains(pdfExportSubDir, 'optResultsOverlap'))
        pdfExportSubDir = strrep(pdfExportSubDir, 'optResultsOverlap', '');
    elseif (contains(pdfExportSubDir, 'optResults'))
        pdfExportSubDir = strrep(pdfExportSubDir, 'optResults', '');
    else
        error('Did not find optResults or optResultsOverlap substring in: ''%s''.', pdfExportSubDir);
    end

    % Generate the path if we need to
    RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
        pdfExportSubDir, theOptimizationResultsFileName, ...
        'generateMissingSubDirs', generateMissingSubDirs);

    pdfExportSubDir = fullfile(pdfExportSubDir, 'optResultsSummaryPDFs');

end
