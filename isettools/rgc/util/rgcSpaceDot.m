function [respCenter, respSurround] = rgcSpaceDot(mosaic, input)
% Spatial inner product of stimulus and each RF.
%
%  [spRespCenter, spRespSurround] = rgcSpaceDot(mosaic,spCenter,spSurround)
%
% Compute the inner product of the center and surround of each receptive
% field with the input.  The input is (typically) a bipolar photocurrent
% image.
%
% This function extracts the relevant x and y coordinates of the stimulus
% at each temporal frame and computes the inner product with the center and
% surround of each cell's RF. The final result is a time series for the
% center and surround of each RGC receptive field.  
%
% The sum of the center and surround is typically the total response,
% though there may be future models in which the center and surround are
% not summed, but rather combined in a more complex formula.
%
% Inputs:
%   mosaic   - rgcMosaic object
%   input    - 3D center input in (x, y, time, color channel) format
%
% Outputs:
%  respCenter:     Response over time from the center
%  respSurround:   Response over time from the surround
%   
% Programming
%   In some cases, the input has several color channels (don't ask) and so
%   we do it for every color channel.  This has to do with the need to
%   handle osDisplayRGB case.
%
% JRG,BW ISETBIO TEAM, 2015

%% init parameters
nSamples = size(input, 3);
nColors  = size(input, 4);
nCells   = mosaic.get('mosaic size');

% pre-allocate space
respCenter   = zeros([nCells(1), nCells(2), nSamples, nColors]);
respSurround = zeros([nCells(1), nCells(2), nSamples, nColors]);

%% Do the inner product.

% The middle cell is at (0,0).  The offset tells us how far offset the 1st
% cell is from (0,0).
switch class(mosaic)
    case 'rgcPhys'
        offset = [0 0];
    otherwise
        % This is the upper leftmost point, I think.
        % offset = [rowConv colConv] .* mosaic.cellLocation{1,1};
        
        micronsToBipolars = mosaic.Parent.col/(1e6*mosaic.Parent.size);
        offset = micronsToBipolars*mosaic.cellLocation{1,1};
end

% BW:  I want to get rid of this nColors element of the loop.
if nColors > 1
    fprintf('nColors bigger than 1.  Call BW and tell him what you are doing\n'); 
end
for cc = 1 : nColors
    for ii = 1 : nCells(1)
        for jj = 1 : nCells(2)
            
            % Get RF of the center and surround of this cell.  These data
            % are not on the mosaic, but rather they are centered at (0,0).
            spRFcenter   = mosaic.sRFcenter{ii, jj};
            spRFsurround = mosaic.sRFsurround{ii, jj};
            % vcNewGraphWin; imagesc(spRFcenter)            
            
            % Center positions of the stimulus used for the inner product
            % with the RF. 
            [stimX, stimY] = mosaic.stimPositions(ii,jj);
            
            % Find the RF location indices, gz, that are within size of
            % stimulus
            gz = find(stimX - offset(1) >= 1 & ...
                stimY - offset(2) >= 1 & ...
                stimX - offset(1) <= size(input,1) & ...
                stimY - offset(2) <= size(input,2) );
            
            % Extract the part of the stimulus that will interact with the
            % RF.  We take this part of the stimulus for all points in
            % time.
            stimV = input(floor(stimX(gz)-offset(1)), ...
                floor(stimY(gz)-offset(2)), :, cc);
            % Visualize the selected portion of the stimulus
            % vcNewGraphWin; ieMovie(stimV);
            
            % Compute and store the inner product of the rf weights and
            % the stimulus computed across all of the time points.
            rfC   = RGB2XWFormat(spRFcenter(gz,gz));
            rfS   = RGB2XWFormat(spRFsurround(gz,gz));
            stimV = RGB2XWFormat(stimV);

            respCenter(ii,jj,:,cc)   = 40*rfC'* stimV;
            respSurround(ii,jj,:,cc) = 40*rfS'* stimV;
            
        end
    end
end

end