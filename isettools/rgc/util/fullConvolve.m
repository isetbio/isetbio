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
% Example:
%    [fullResponse, nlResponse] = fullConvolve(rgc1.mosaic{1}, spResponseCenter, spResponseSurround);
%     
% (c) isetbio
% 09/2015 JRG

% Find bounds for size of input and output
spResponseSize = size(spResponseCenter{1,1}(:,:,1,1));
nSamples = size(spResponseCenter{1,1},3);
channelSize = size(spResponseCenter{1,1},4);

nCells = size(mosaic.cellLocation);

fprintf('Temporal Convolution, %s:     \n', mosaic.cellType);

% set flag to keep track of whether each cell has its own temporal impulse
% response function
if length(mosaic.tCenter) == 3
    lenflag = 0;
else 
    lenflag = 1;
end


tic
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        for rgbIndex = 1:channelSize
            
            % fullResponseRSRGB = zeros(size(spResponseCenter{1,1}));
            
            % Get temporal impulse response functions
            if lenflag == 0
                temporalIRCenter = mosaic.tCenter{rgbIndex};
                temporalIRSurround = mosaic.tSurround{rgbIndex};
            else
                
                temporalIRCenter = mosaic.tCenter{xcell,ycell};
                temporalIRSurround = mosaic.tSurround{xcell,ycell};
            end
            
            % Reshape the spatial responses from spConvolve to allow for
            % efficient computation of the convolution with the temp IRF
            
            spResponseSize = size(spResponseCenter{xcell,ycell}(:,:,:,rgbIndex));            
            
            spResponseCenterRS = reshape(squeeze(spResponseCenter{xcell,ycell}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
            spResponseSurroundRS = reshape(squeeze(spResponseSurround{xcell,ycell}(:,:,:,rgbIndex)), spResponseSize(1)*spResponseSize(2), nSamples);
    
            if (sum(temporalIRCenter(:)-temporalIRSurround(:)) == 0) 
                % if the temporal impulse responses for center and surround are the same, combine before convolution for efficiency                                             
                fullResponseRSCombined = convn(spResponseCenterRS-spResponseSurroundRS, temporalIRCenter','full');
                
                % Specify starting and ending time coordinates
                startPoint = 1; endPoint = nSamples+length(temporalIRCenter)-1;
                fullResponseRSRGB = zeros(spResponseSize(1)*spResponseSize(2),length(startPoint:endPoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCombined(:,startPoint:endPoint);
                                
            else            
                % if the temporal impulse responses for center and surround are different, convolve both               
                fullResponseRSCenter = convn(spResponseCenterRS, temporalIRCenter','full');
                fullResponseRSSurround = convn(spResponseSurroundRS, temporalIRSurround','full');
                
                % Specify starting and ending time coordinates
%                 startPoint = length(temporalIRCenter)-1; endPoint = nSamples+length(temporalIRCenter)-1;
                startPoint = 1; endPoint = nSamples+length(temporalIRCenter)-1;
                % Take difference between center and surround response                
                fullResponseRSRGB = zeros(spResponseSize(1)*spResponseSize(2),length(startPoint:endPoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,startPoint:endPoint) - fullResponseRSSurround(:,startPoint:endPoint);
                
            end            
        end      
                        
        % Take the mean of the spatial response over (x,y) at a particular
        % time frame t for each cell
        if isa(mosaic,'rgcMosaicSubunit');
            % For the subunit model, apply the nonlinearity before taking the mean
            genFunction = mosaicGet(mosaic, 'generatorFunction');
            fullResponseRS = sum(genFunction(fullResponseRSRGB),3);            
            fullResponse{xcell,ycell,1} = mean(fullResponseRS);            
            nlResponse{xcell,ycell} = (mean(fullResponseRS,1));
        elseif isa(mosaic,'rgcMosaicPhys')
            
            % For all other models, apply the nonlinearity after
            fullResponseRS = sum(fullResponseRSRGB,3);                     
            fullResponse{xcell,ycell,1} = sum(fullResponseRS) + mosaic.tonicDrive{xcell,ycell};       % mean? sum in ej's code
            % % fullResponse for RGB
            fullResponse{xcell,ycell,2} =  reshape(fullResponseRSRGB, spResponseSize(1), spResponseSize(2), size(fullResponseRSRGB,2), size(fullResponseRSRGB,3));
            genFunction = mosaicGet(mosaic, 'generatorFunction');
            nlResponse{xcell,ycell} = genFunction(sum(fullResponseRS,1) + mosaic.tonicDrive{xcell,ycell});
            
        else
            % For all other models, apply the nonlinearity after
            fullResponseRS = sum(fullResponseRSRGB,3);                     
            fullResponse{xcell,ycell,1} = mean(fullResponseRS);       % this is the only difference from the elseif block
            % % fullResponse for RGB
            fullResponse{xcell,ycell,2} =  reshape(fullResponseRSRGB, spResponseSize(1), spResponseSize(2), size(fullResponseRSRGB,2), size(fullResponseRSRGB,3));
            if ~isa(mosaic, 'rgcMosaicLinear'); 
                genFunction = mosaicGet(mosaic, 'generatorFunction');
                nlResponse{xcell,ycell} = genFunction(mean(fullResponseRS,1));
            else
                
                nlResponse{xcell,ycell} = 0;
            end
        end
        
       
                                      
    end
end
toc