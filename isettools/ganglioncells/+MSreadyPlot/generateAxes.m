function theAxes = generateAxes(hFig,ff, varargin)

    p = inputParser;
    p.addParameter('figPosition', [], @(x)(isempty(x)||(numel(x) == 2)));
    p.parse(varargin{:});
    figPosition = p.Results.figPosition;

    if (isempty(figPosition))
        set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
    else
        set(hFig, 'Position', [figPosition(1) figPosition(2) ff.figureSize(1) ff.figureSize(2)], 'Color', [1 1 1]);
    end

    % Generate axes
    for iRow = 1:size(ff.subplotPosVectors,1)
        for iCol = 1:size(ff.subplotPosVectors,2)
            figure(hFig);
            theAxes{iRow,iCol} = subplot('Position', ff.subplotPosVectors(iRow, iCol).v);
        end
    end

end
