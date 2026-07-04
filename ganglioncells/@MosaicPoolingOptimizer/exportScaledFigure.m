function exportScaledFigure(sourcePDF, destinationPDF, scaleFactor)

    cpdfPath = '/Users/nicolas/Documents/Manuscripts';
    cpdfPath = fullfile(cpdfPath,'CPDF/cpdf');

    if (~exist(cpdfPath, 'file'))
        error('cpdf software not found in expected location: %s', cpdfPath);
    end

    sourcePDF = strrep(sourcePDF, ' ', '\ ');
    destinationPDF = strrep(destinationPDF, ' ', '\ ');

    fprintf('Scaled version of %s saved to %s\n', sourcePDF, destinationPDF);
    commandString = sprintf('%s -scale-page "%2.2f %2.2f" %s -o %s', ...
        cpdfPath, scaleFactor, scaleFactor, sourcePDF, destinationPDF);
    system(commandString);
end