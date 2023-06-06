function theAxes = generateAxes(hFig,ff)

    set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
    
    % Generate axes
    for iRow = 1:size(ff.subplotPosVectors,1)
        for iCol = 1:size(ff.subplotPosVectors,2)
            figure(hFig);
            theAxes{iRow,iCol} = subplot('Position', ff.subplotPosVectors(iRow, iCol).v);
        end
    end

end
