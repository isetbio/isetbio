function exportScaledFigure(sourcePDF, destinationPDF, scaleFactor)

    cpdfCommandPath = '/Users/nicolas/Documents/Developer/Matlab/CPDF/cpdf';
    sourcePDF = strrep(sourcePDF, ' ', '\ ');
    destinationPDF = strrep(destinationPDF, ' ', '\ ');

    fprintf('Scaled version of %s saved to %s\n', sourcePDF, destinationPDF);
    commandString = sprintf('%s -scale-page "%2.2f %2.2f" %s -o %s', ...
        cpdfCommandPath, scaleFactor, scaleFactor, sourcePDF, destinationPDF);
    system(commandString);
end