function [respCenter, respSurround] = rgcSpaceDot(mosaic, input)
% Spatial inner product of stimulus and each RF.
%
%  [spRespCenter, spRespSurround] = rgcSpaceDot(mosaic,input)
%
% Compute the inner product of the center and surround of each receptive
% field with the input.  The input is (typically) a bipolar mosaic with a
% photocurrent response.
%
% This function extracts the x and y coordinates of the input at each
% temporal frame and computes the inner product with the center and
% surround of each cell's RF. The final result is a time series for the
% center and surround for each RGC receptive field.
%
% The sum of the center and surround is the total response, though there
% may be future models in which the center and surround are not summed, but
% rather combined in a more complex formulae.
%
% The cell locations of the RGC are in coordinates of microns with (0,0) at
% the center.  We need to co-register the inputs (which are bipolar cells)
% with the RGC receptive fields.  The inputs just come in as a 3D matrix,
% without any coordinates. The first stage of the calculation assigns a
% coordinate to each bipolar.
%
% Inputs:
%   mosaic   - rgcMosaic object
%   input    - 3D input in (x, y, time) format
%
% Outputs:
%  respCenter:     Response over time from the center
%  respSurround:   Response over time from the surround
%
% See below for PROGRAMMING TODO
%
% JRG,BW ISETBIO TEAM, 2015

%% init parameters
nSamples = size(input, 3);
nCells   = mosaic.get('mosaic size');

% pre-allocate space
respCenter   = zeros([nCells(1), nCells(2), nSamples]);
respSurround = zeros([nCells(1), nCells(2), nSamples]);

%% Do the inner product.

% The middle cell is at (0,0).  The offset tells us how far offset the 1st
% cell is from (0,0). This value is the same for all cells in the mosaic so
% computed outside of the loop.
patchSizeUM = 1e6*mosaic.Parent.size;   % In microns
bipolarsPerMicron = size(input(:,:,1)) ./ patchSizeUM;   % cells/micron
nTime = size(input,3);

for ii = 1 : nCells(1)
    for jj = 1 : nCells(2)
        
        % Get RF of the center and surround of this cell.  These data
        % are not on the mosaic, but rather they are centered at (0,0).
        spRFcenter   = mosaic.sRFcenter{ii, jj};
        spRFsurround = mosaic.sRFsurround{ii, jj};
        % vcNewGraphWin; imagesc(spRFcenter)
        
        % Row and col positions of the input used for the inner product
        % with the RGC RF. The row/col values are sample positions of the
        % bipolar input. 
        [inputRow, inputCol] = mosaic.inputPositions(ii,jj,bipolarsPerMicron);
        
        % Find the rows and columns within the stimulus range
        nRow = length(inputRow); nCol = length(inputCol);
        mxRow = size(input,1); mxCol = size(input,2);
        inputValues = zeros(nRow,nCol,nTime);
        
        % Sometimes the RF extends outside of the size of the input.  So we
        % clip it here to keep the values within the input range, these are
        % the good rows and cols.
        goodRowInds = (inputRow > 0) & (inputRow < mxRow);
        goodRows = inputRow(goodRowInds);
        goodColInds = (inputCol > 0) & (inputCol < mxCol);
        goodCols = inputCol(goodColInds);
        inputValues(goodRowInds,goodColInds,:) = input(goodRows,goodCols,:);
        
        % The inputValues are a time series.  So we put it into XW (which
        % means space-time) formatting.
        inputValues   = RGB2XWFormat(inputValues);
        rfC           = spRFcenter(:);
        rfS           = spRFsurround(:);
        
        % Put the RF center as a row and multiply it by the XT matrix
        respCenter(ii,jj,:)   = rfC' * inputValues;
        respSurround(ii,jj,:) = rfS' * inputValues;
        
        % PROGRAMMING TODO
        % We have to handle edge cases.  Let's see if we can improve on
        % this particular 'find' function. For example, ...
        % 
        % gz = (inputRow > 1) && (inputRow < mx)
        %
        % Find the RF location indices, gz, that are within the size of
        % stimulus
        %         gz = find(inputRow - offset(1) >= 1 & ...
        %             inputCol - offset(2) >= 1 & ...
        %             inputRow - offset(1) <= size(input,1) & ...
        %             inputCol - offset(2) <= size(input,2) );
        %
        %         % Extract the part of the input that will interact with the
        %         % RF.  We take this part of the input for all points in
        %         % time.
        %         inputValues = input(floor(inputRow(gz)-offset(1)), ...
        %             floor(inputCol(gz)-offset(2)), :);
        %         % Visualize the selected portion of the stimulus
        %         % vcNewGraphWin; ieMovie(stimV);
        %
        %         % Compute and store the inner product of the rf weights and
        %         % the stimulus computed across all of the time points.
        %         rfC   = RGB2XWFormat(spRFcenter(gz,gz));
        %         rfS   = RGB2XWFormat(spRFsurround(gz,gz));
        %         inputValues = RGB2XWFormat(inputValues);
        %
        %         % What is this 40 doing here? PROGRAMMING TODO
        %         respCenter(ii,jj,:)   = 40*rfC'* inputValues;
        %         respSurround(ii,jj,:) = 40*rfS'* inputValues;
        
    end
end


end