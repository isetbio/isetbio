function [respCenter, respSurround] = rgcSpaceDot(rgcM, input)
% Spatial inner product of stimulus and each RF.
%
% Syntax:
%   [spRespCenter, spRespSurround] = rgcSpaceDot(mosaic, input)
%
%    Compute the inner product of the center and surround of each receptive
%    field with the input. The input is (typically) a bipolar mosaic with a
%    photocurrent response.
%
%    This function extracts the x and y coordinates of the input at each
%    temporal frame and computes the inner product with the center and
%    surround of each cell's RF. The final result is a time series for the
%    center and surround for each RGC receptive field.
%
%    The sum of the center and surround is the total response, though there
%    may be future models in which the center and surround are not summed,
%    but rather combined in a more complex formulae.
%
%    The cell locations of the RGC are in coordinates of microns with (0,
%    0) at the center. We need to co-register the inputs (which are bipolar
%    cells) with the RGC receptive fields. The inputs just come in as a 3D
%    matrix, without any coordinates. The first stage of the calculation
%    assigns a coordinate to each bipolar.
%
% Inputs:
%    mosaic       - Object. A rgcMosaic object
%    input        - Matrix. A 3D input matrix in (x, y, time) format.
%
% Outputs:
%    respCenter   - Matrix. A 3D matrix comprised of the response over time
%                   from the center.
%    respSurround - Matric. A 3D matrix comprised of the response over time
%                   from the surround.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/15  JRG, BW  ISETBIO TEAM, 2015
%    05/31/19  JNM      Documentation pass

%% init parameters
nSamples = size(input, 3);            % Temporal samples
nCells = rgcM.get('mosaic samples');  % RGC cells

% pre-allocate space
respCenter = zeros([nCells(1), nCells(2), nSamples]);
respSurround = zeros([nCells(1), nCells(2), nSamples]);

%% Do the inner product.
% The middle cell is at (0, 0). The offset tells us how far offset the 1st
% cell is from (0, 0). This value is the same for all cells in the mosaic
% so computed outside of the loop.
% patchSizeUM = 1e6 * mosaic.parent.size;                   % In microns
% bipolarsPerMicron = size(input(:, :, 1)) ./ patchSizeUM;  % cells/micron
nTime = size(input, 3);

for ii = 1 : nCells(1)
    for jj = 1 : nCells(2)
        % Get RF of the center and surround of this cell. These data
        % are not on the mosaic, but rather they are centered at (0, 0).
        spRFcenter = rgcM.sRFcenter{ii, jj};
        spRFsurround = rgcM.sRFsurround{ii, jj};
        % vcNewGraphWin; imagesc(spRFcenter)

        % Row and col positions of the input used for the inner product
        % with the RGC RF. The row/col values are sample positions of the
        % bipolar input. We need a unit test to show that these positions
        % change properly as we step through the RGC positions.
        [inputRow, inputCol] = rgcM.inputPositions(ii, jj);

        % Find the rows and columns within the stimulus range
        nRow = length(inputRow); nCol = length(inputCol);
        mxRow = size(input, 1);    mxCol = size(input, 2);
        inputValues = zeros(nRow, nCol, nTime);

        % Sometimes the RF extends outside of the size of the input. So we
        % clip it here to keep the values within the input range, these are
        % the good rows and cols.
        goodRowInds = (inputRow > 0) & (inputRow < mxRow);
        goodRows = inputRow(goodRowInds);
        goodColInds = (inputCol > 0) & (inputCol < mxCol);
        goodCols = inputCol(goodColInds);
        inputValues(goodRowInds, goodColInds, :) = ...
            input(goodRows, goodCols, :);

        % badRowInds = (inputRow < 0) | (inputRow > mxRow);
        % badRows = inputRow(badRowInds);
        % badColInds = (inputCol < 0) | (inputCol > mxCol);
        % badCols = inputCol(badColInds);

        % inputGood = RGB2XWFormat(input(goodRows, goodCols, :));
        % meanInputGood(1, 1, :) = mean(inputGood);
        % inputBad = repmat(meanInputGood, ...
        %     [length(badRows), length(badCols), 1]);
        % inputValues(badRowInds, badColInds, :) = ...
        %     inputBad;%input(badRows, badCols, :);

        % The inputValues are a time series. So we put it into XW (which
        % means space-time) formatting.
        inputValues = RGB2XWFormat(inputValues);
        rfC = spRFcenter(:);
        rfS = spRFsurround(:);

        % Put the RF center as a row and multiply it by the XT matrix
        respCenter(ii, jj, :) = rfC' * inputValues;
        respSurround(ii, jj, :) = rfS' * inputValues;
    end
end

end