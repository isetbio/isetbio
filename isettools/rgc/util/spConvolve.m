function [spRespCenter, spRespSurround] = spConvolve(mosaic, spCenter, spSurround)
% This is a util function of the @rgc parent class, for a separable
% STRF finds the 2D convolution of the spatial RF for every time frame.
% 
%  [spRespCenter, spRespSurround] = spConvolve(mosaic,spCenter,spSurround)
%   
% This function extracts the relevant x and y coordinates of each temporal
% frame of the stimulus and convolves that 2D image with the spatial center
% and surround RFs of each cell for all temporal frames.
% 
% Inputs:
%   mosaic     - rgcMosaic object
%   spCenter   - 3D center input in (x, y, t) format
%   spSurround - 3D surround input in (x, y, t) format
% 
% Outputs: 
%   spatial response of the 2D center and surround receptive fields.
%
%
% JRG, ISETBIO TEAM, 2015

% init parameters
nSamples = size(spCenter, 3);
channelSize = size(spCenter, 4);
nCells = size(mosaic.cellLocation);

% pre-allocate space
spRespCenter = cell(nCells);
spRespSurround = cell(nCells);

for rgbIndex = 1 : channelSize
    for ii = 1 : nCells(1)
        for jj = 1 : nCells(2)
            % Get RF 2D images
            spRFcenter = mosaic.sRFcenter{ii, jj};
            spRFsurround = mosaic.sRFsurround{ii, jj};
            
            spRespCenter{ii,jj} = zeros([size(spRFcenter) nSamples channelSize]);
            spRespSurround{ii,jj} = zeros([size(spRFsurround) nSamples channelSize]);
            
            % Get cell location
            stimCenterCoords = mosaic.cellLocation{ii,jj};
            switch class(mosaic)
                case 'rgcPhys'
                    % Find the spatial extent of the RF in numbers of rfDiameter
                    extent = round(size(mosaic.sRFcenter{1,1},1)/mosaic.rfDiameter);
                    offset = [0 0];
                
                    % Should rfDimater be size RF center? Yes!
                    stimX =  ceil((stimCenterCoords(1) - floor(extent/2*mosaic.rfDiameter))):floor((stimCenterCoords(1) + floor((extent/2)*mosaic.rfDiameter )));
                    stimY =  ceil((stimCenterCoords(2) - floor(extent/2*mosaic.rfDiameter))):floor((stimCenterCoords(2) + floor((extent/2)*mosaic.rfDiameter )));
                case 'rgcSubunit'
                    extent = 1;
                    if mosaic.cellLocation{1,1}(1) > 0
                        offset(1) = ceil(mosaic.cellLocation{1,1}(1));
                    else
                        offset(1) = floor(mosaic.cellLocation{1,1}(1));
                    end
                
                    if mosaic.cellLocation{1,1}(2) > 0
                        offset(2) = ceil(mosaic.cellLocation{1,1}(2));
                    else
                        offset(2) = floor(mosaic.cellLocation{1,1}(2));
                    end
                    
                    % Should rfDimater be size RF center? Yes!
                    stimX =  floor((stimCenterCoords(1) - floor(extent/2*size(mosaic.sRFcenter{1,1},1)))) : floor((stimCenterCoords(1) + floor((extent/2)*size(mosaic.sRFcenter{1,1},1) )));
                    stimY =  floor((stimCenterCoords(2) - floor(extent/2*size(mosaic.sRFcenter{1,1},2)))) : floor((stimCenterCoords(2) + floor((extent/2)*size(mosaic.sRFcenter{1,1},2))));
                    
                otherwise
                    extent = 1;
                    rowConv = 1; colConv = 1;
                    offset = [rowConv colConv].*mosaic.cellLocation{1,1};
                    
                    % Should rfDimater be size RF center? Yes!
                    stimX =  ceil((stimCenterCoords(1) - floor(extent/2*size(mosaic.sRFcenter{1,1},1)))) : floor((stimCenterCoords(1) + floor((extent/2)*size(mosaic.sRFcenter{1,1},1) )));
                    stimY =  ceil((stimCenterCoords(2) - floor(extent/2*size(mosaic.sRFcenter{1,1},2)))) : floor((stimCenterCoords(2) + floor((extent/2)*size(mosaic.sRFcenter{1,1},2))));
                    
            end
            
            % Ensure indices are within size of stimulus
            gz = find(stimX - offset(1) >= 1 & ...
                stimY - offset(2) >= 1 & ...
                stimX - offset(1) <= size(spCenter,1) & ...
                stimY - offset(2) <= size(spCenter,2) );
            
            spStimCenterV = spCenter(floor(stimX(gz)-offset(1)), ...
                floor(stimY(gz)-offset(2)), :,rgbIndex);
            spStimSurrV = spSurround(floor(stimX(gz)-offset(1)), ...
                floor(stimY(gz)-offset(2)), :, rgbIndex);
            
            switch class(mosaic)
                case 'rgcSubunit'
                    spRC = bsxfun(@times,spRFcenter(gz,gz),spStimCenterV);
                    spRS = bsxfun(@times,spRFsurround(gz,gz),spStimCenterV);
                    spRespCenter{ii,jj}(gz, gz, :, rgbIndex) = ...
                        10 * mosaic.rectifyFunction(spRC) / length(gz)^2;
                    spRespSurround{ii,jj}(gz,gz, :,rgbIndex) = ...
                        10 * mosaic.rectifyFunction(spRS) / length(gz)^2;
                case 'rgcPhys'
                    spRespCenter{ii,jj}(gz, gz, :, rgbIndex) = ...
                        bsxfun(@times, spRFcenter(gz, gz), ...
                        spStimCenterV-spStimSurrV);
                otherwise  % other types of rgc, include LNP, etc.
                    spRespCenter{ii,jj}(gz, gz, :, rgbIndex) = ...
                        bsxfun(@times, spRFcenter(gz,gz), spStimCenterV);
                    spRespSurround{ii,jj}(gz, gz, :, rgbIndex) = ...
                        bsxfun(@times, spRFsurround(gz,gz), spStimCenterV);
            end
        end
    end
end

end