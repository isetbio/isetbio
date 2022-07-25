function appendFigAsSeparatePageInPDFdoc(pdfFileName,figHandle,dpi)

        % If no handle is provided, use the current figure as default
    if nargin<1
        [fileName,pathName] = uiputfile('*.pdf','Save to PDF file:');
        if fileName == 0; return; end
        pdfFileName = [pathName,fileName];
    end
    if nargin<2
        handle = gcf;
    end
    if nargin<3
        dpi = 150;
    end
        
    % Backup previous settings
    prePaperType = get(figHandle,'PaperType');
    prePaperUnits = get(figHandle,'PaperUnits');
    preUnits = get(figHandle,'Units');
    prePaperPosition = get(figHandle,'PaperPosition');
    prePaperSize = get(figHandle,'PaperSize');

    % Make changing paper type possible
    set(figHandle,'PaperType','<custom>');

    % Set units to all be the same
    set(figHandle,'PaperUnits','inches');
    set(figHandle,'Units','inches');

    % Set the page size and position to match the figure's dimensions
    paperPosition = get(figHandle,'PaperPosition');
    position = get(figHandle,'Position');
    set(figHandle,'PaperPosition',[0,0,position(3:4)]);
    set(figHandle,'PaperSize',position(3:4));

    set(figHandle,'InvertHardCopy','off');
    
    if (exist(pdfFileName, 'file'))
        tmpFileName = sprintf('tmpPDF.pdf');
        mergedFileName = sprintf('mergedPDF.pdf');
        print(figHandle,'-dpdf', tmpFileName,sprintf('-r%d',dpi));
        status = system(sprintf('/opt/homebrew/bin/pdfunite %s %s %s', pdfFileName, tmpFileName, mergedFileName));
        if (status ~= 0)
            fprintf('Error during pdfunite. Have you installed poppler (''brew install poppler'')?');
            system(sprintf('rm %s', tmpFileName));
        else
            system(sprintf('rm %s', tmpFileName));
            system(sprintf('mv %s %s', mergedFileName, pdfFileName));
            fprintf('\nNicePlot: PDF appended as the last page of %s.\n', pdfFileName);
        end
    else
        fprintf('\nNicePlot: PDF saved as the first page of %s.\n', pdfFileName);
        print(figHandle,'-dpdf', pdfFileName,sprintf('-r%d',dpi));
    end

    % Restore the previous settings
    set(figHandle,'PaperType',prePaperType);
    set(figHandle,'PaperUnits',prePaperUnits);
    set(figHandle,'Units',preUnits);
    set(figHandle,'PaperPosition',prePaperPosition);
    set(figHandle,'PaperSize',prePaperSize);
    
end

