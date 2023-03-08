function theAxes = generateAxes(hFig,ff)

% Generate axes
for iRow = 1:size(ff.subplotPosVectors,1)
    for iCol = 1:size(ff.subplotPosVectors,2)
        figure(hFig);
        theAxes{iRow,iCol} = subplot('Position', ff.subplotPosVectors(iRow, iCol).v);
    end
end

end
