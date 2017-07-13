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
% JRG,BW ISETBIO TEAM, 2015

%% init parameters
nSamples = size(input, 3);
nCells   = mosaic.get('mosaic size');

% pre-allocate space
respCenter   = zeros([nCells(1), nCells(2), nSamples]);
respSurround = zeros([nCells(1), nCells(2), nSamples]);

%% Do the inner product.

% The middle cell is at (0,0).  The offset tells us how far offset the 1st
% cell is from (0,0).
switch class(mosaic)
    case 'rgcPhys'
        offset = [0 0];
    otherwise
        % BW is confused by this calculation.  The cellLocation term is in
        % samples on the cone mosaic.  What is being computed here appears
        % to be the position in microns of the bipolar samples.  Is that
        % right?
        
        % This is the upper leftmost point, I think.
        % offset = [rowConv colConv] .* mosaic.cellLocation{1,1};
        
        % I think this section should be replaced by
        %
        %    1 / bp.mosaic{1}.rfSize('units','um');
        %
        % The same code appears elswhere (rgcMosaic.stimPositions)
        % (BW)
        
        % This code replaces something JRG had that confused me.  That code
        % as (size/(1e6*size)), which is always 1e6.  I think what is
        % intended is the number of bipolar receptive fields per micron.
        % When the RF is 2 microns, this is 1/2.  (BW).
        patchSizeUM = 1e6*mosaic.Parent.size;   % In microns
        
        % The bipolar mosaics at this point are all the same row/col count.
        % But they may not be in the future.  So, what do we do about that?
        bpRowCol = size(input);
        
        bipolarsPerMicron = bpRowCol(1:2) ./ patchSizeUM;   % cells/micron
        
        % The cellLocation positions are in microns and centered at 0,0.
        % The offset is the
        offset = bipolarsPerMicron .* mosaic.cellLocation{1,1};
end

for ii = 1 : nCells(1)
    for jj = 1 : nCells(2)
        
        % Get RF of the center and surround of this cell.  These data
        % are not on the mosaic, but rather they are centered at (0,0).
        spRFcenter   = mosaic.sRFcenter{ii, jj};
        spRFsurround = mosaic.sRFsurround{ii, jj};
        % vcNewGraphWin; imagesc(spRFcenter)
        
        % Center positions of the stimulus used for the inner product
        % with the RF. stimPositions are indices of the bipolar input. We
        % might rename this to be inputPositions
        [inputX, inputY] = mosaic.stimPositions(ii,jj,bipolarsPerMicron);
        
        % Find the RF location indices, gz, that are within size of
        % stimulus
        gz = find(inputX - offset(1) >= 1 & ...
            inputY - offset(2) >= 1 & ...
            inputX - offset(1) <= size(input,1) & ...
            inputY - offset(2) <= size(input,2) );
        
        % Extract the part of the input that will interact with the
        % RF.  We take this part of the input for all points in
        % time.
        stimV = input(floor(inputX(gz)-offset(1)), ...
            floor(inputY(gz)-offset(2)), :);
        % Visualize the selected portion of the stimulus
        % vcNewGraphWin; ieMovie(stimV);
        
        % Compute and store the inner product of the rf weights and
        % the stimulus computed across all of the time points.
        rfC   = RGB2XWFormat(spRFcenter(gz,gz));
        rfS   = RGB2XWFormat(spRFsurround(gz,gz));
        stimV = RGB2XWFormat(stimV);
        
        respCenter(ii,jj,:)   = 40*rfC'* stimV;
        respSurround(ii,jj,:) = 40*rfS'* stimV;
        
    end
end


end