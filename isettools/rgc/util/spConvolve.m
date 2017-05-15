function [respCenter, respSurround] = spConvolve(mosaic, input)
% A separable space-time 2D convolution for center/surround
%
%  [spRespCenter, spRespSurround] = spConvolve(mosaic,spCenter,spSurround)
%
% Compute the inner product of the center and surround of each receptive
% field with the input, which is typically a bipolar photocurrent image.
%
% This function extracts the relevant x and y coordinates of the stimulus
% at each temporal frame and computes the inner product with the center and
% surround components of the RF. The final result is a time series for the
% center and surround of each RGC receptive field.  The sum of these is
% typically considered the total response, though there may be future
% models in which the center and surround are not summed.
%
% Inputs:
%   mosaic   - rgcMosaic object
%   input    - 3D center input in (x, y, time, color channel) format
%
% Outputs:
%  respCenter:     Response over time from the center
%  respSurround:   Response over time from the surround
%   
%
%   N.B. In some cases, the input has several color channels (don't ask) and so
%   we do it for every color channel.  This has to do with the need to
%   handle osDisplayRGB case.
%
% JRG, ISETBIO TEAM, 2015
% 7/2016, JRG/BW updated

%% init parameters
nSamples = size(input, 3);
nColors  = size(input, 4);
nCells   = mosaic.get('mosaic size');

% pre-allocate space
% @JRG - Could turn this into a matrix, not a cell array
respCenter   = zeros([nCells(1), nCells(2), nSamples, nColors]);
respSurround = zeros([nCells(1), nCells(2), nSamples, nColors]);

%% Do the convolution.

% BW:  I don't understand this. Am deleting because it doesn't do anything
% anyway.
% rowConv = 1; colConv = 1;

% The middle cell is at (0,0).  This tells us how far offset the 1st cell
% is from (0,0).
switch class(mosaic)
    case 'rgcPhys'
        offset = [0 0];
    otherwise
        % offset = [rowConv colConv].*mosaic.cellLocation{1,1};
        offset = mosaic.cellLocation{1,1};
end

% BW:  I want to get rid of this nColors thing.
for cc = 1 : nColors
    for ii = 1 : nCells(1)
        for jj = 1 : nCells(2)
            % Get RF of the center and surround of this cell
            spRFcenter   = mosaic.sRFcenter{ii, jj};
            spRFsurround = mosaic.sRFsurround{ii, jj};
            % vcNewGraphWin; imagesc(spRFcenter)            
            
            % We want a routine that pulls out the relevant portion of the
            % input for this RF.  This code is a bit hard to read,
            % so we should simplify.
            
            % Positions of the stimulus used for the inner product
            [stimX, stimY] = mosaic.stimPositions(ii,jj);
            
            % Find the RF location indices that are within size of stimulus
            gz = find(stimX - offset(1) >= 1 & ...
                stimY - offset(2) >= 1 & ...
                stimX - offset(1) <= size(input,1) & ...
                stimY - offset(2) <= size(input,2) );
            
            % Clip the proper part of the stimulus with respect to the RF
            % We do all of the time dimension in a single matrix
            stimV = input(floor(stimX(gz)-offset(1)), ...
                floor(stimY(gz)-offset(2)), :, cc);
            % Visualize the selected portion of the stimulus
            % vcNewGraphWin; ieMovie(stimV);
            
            % This might be quick, or perhaps this calls for a bsxfun()
            % routine to do the multiply.
            
            % Apply the rf weights to the stimulus across all of the time
            % points.
            rfC = RGB2XWFormat(spRFcenter(gz,gz));
            rfS = RGB2XWFormat(spRFsurround(gz,gz));
            stimV = RGB2XWFormat(stimV);
            
            respCenter(ii,jj,:,cc)   = rfC'* stimV;
            respSurround(ii,jj,:,cc) = rfS'* stimV;
            
        end
    end
end

end