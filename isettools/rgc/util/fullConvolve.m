function [fullResponse, nlResponse] = fullConvolve(mosaic, spResponseCenter, spResponseSurround)
% fullConvolve: a util function of the @rgc parent class, for a separable
% STRF finds the 1D convolution of the temporal impulse response with the
% output signal of the spatial convolution operation.
% 
%         [fullResponse, nlResponse] = fullConvolve(mosaic, spResponseCenter, spResponseSurround);
%     
% Inputs: the mosaic object, the center and surround spatial resposes for
%   each temporal stimulus frame, from spConvolve.m.
% 
% Outputs: the full response for each cell over t frames and the nonlinear
%   response following the generator lookup function.
%     
% JRG, ISETBIO TEAM, 2015

% Find bounds for size of input and output
nSamples = size(spResponseCenter{1,1},3);
channelSize = size(spResponseCenter{1,1},4);

nCells = size(mosaic.cellLocation);

% set flag to keep track of whether each cell has its own temporal impulse
% response function
if length(mosaic.tCenter) == 3
    lenflag = 0;
else 
    lenflag = 1;
end

fullResponse = cell([nCells 2]);
nlResponse = cell(nCells);

for ii = 1:nCells(1)
    for jj = 1:nCells(2)
        for rgbIndex = 1:channelSize
            % Get temporal impulse response functions
            if lenflag == 0
                temporalIRCenter = mosaic.tCenter{rgbIndex};
                temporalIRSurround = mosaic.tSurround{rgbIndex};
            else
                temporalIRCenter = mosaic.tCenter{ii, jj};
                temporalIRSurround = mosaic.tSurround{ii, jj};
            end
            
            % Reshape the spatial responses from spConvolve to allow for
            % efficient computation of the convolution with the temp IRF
            [spCenter,r,c] = RGB2XWFormat(spResponseCenter{ii,jj}(:,:,:,rgbIndex));
            spSurround = RGB2XWFormat(spResponseSurround{ii,jj}(:,:,:,rgbIndex));
    
            if (sum(temporalIRCenter(:)-temporalIRSurround(:)) == 0) 
                % if the temporal impulse responses for center and surround
                % are the same, combine before convolution for efficiency
                
                % Convolution implementation
                fullResponseRSCombined = convn(spResponseCenterRS-spResponseSurroundRS, temporalIRCenter','full');
                
                if isa(mosaic,'rgcSubunit');
                    scaleFactor = 73.73;
                else
                    scaleFactor = 1;
                end
                
                % make the response zero mean
                fullResponseRSCombined = scaleFactor * bsxfun(@minus, ...
                    fullResponseRSCombined,mean(fullResponseRSCombined,2));
               
                % Specify starting and ending time coordinates
                sPoint = 1; ePoint = nSamples;
                fullResponseRSRGB = zeros(r*c,length(sPoint:ePoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCombined(:,sPoint:ePoint);
                                
            else
                fullCenter = convn(spCenter, temporalIRCenter','full');
                fullSurround = convn(spSurround, temporalIRSurround','full');
                
                % Specify starting and ending time coordinates
                sPoint = 1; ePoint = nSamples;
                
                % Take difference between center and surround response                
                fullResponseRSRGB = zeros(r*c,length(sPoint:ePoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullCenter(:,sPoint:ePoint) - fullSurround(:,sPoint:ePoint);
            end            
        end      
                        
        if isa(mosaic,'rgcPhys')
            % For all other models, apply the nonlinearity after
            fullResponseRS = sum(fullResponseRSRGB,3);                     
            fullResponse{ii,jj,1} = sum(fullResponseRS) + mosaic.tonicDrive{ii,jj};       % mean? sum in ej's code
            % % fullResponse for RGB
            fullResponse{ii,jj,2} =  reshape(fullResponseRSRGB, r, c, size(fullResponseRSRGB,2), []);
            
        else
            % For all other models, apply the nonlinearity after
            fullResponseRS = sum(fullResponseRSRGB, 3);
            
            % Need to normalize by 13x13 RF size in lab physio code
            fullResponse{ii,jj,1} = sum(fullResponseRS) + mosaic.tonicDrive{ii,jj};   
            
            % fullResponse for RGB
            fullResponse{ii,jj,2} =  reshape(fullResponseRSRGB, r, c, size(fullResponseRSRGB,2), []);
            if ~isa(mosaic, 'rgcLinear'); 
                genFunction = mosaicGet(mosaic, 'generatorFunction');
                nlResponse{ii,jj} = genFunction(mean(fullResponseRS,1));
            else
                nlResponse{ii,jj} = 0;
            end
        end                      
    end
end

end