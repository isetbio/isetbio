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
                
                % This line was used before the fft was implemented
%                 fullResponseRSCombined = convn(spResponseCenterRS-spResponseSurroundRS, temporalIRCenter','full');
% %                 fullResponseRSCombined = ifft(fft(spResponseCenterRS-spResponseSurroundRS).*fft(temporalIRCenter'));

%%%% FFT that works
                spResponseCenterRSp = [spResponseCenterRS];% zeros([size(spResponseCenterRS,1) size(temporalIRCenter,1)])];
                spResponseSurroundRSp = [spResponseSurroundRS];% zeros([size(spResponseSurroundRS,1) size(temporalIRCenter,1)])];
                temporalIRCenterp = repmat([temporalIRCenter' zeros(1,[-size(temporalIRCenter,1)+size(spResponseCenterRS,2)])],size(spResponseCenterRS,1) ,1);
                
                fullResponseRSCombined = ifft(fft(spResponseCenterRSp'-spResponseSurroundRSp').*fft(temporalIRCenterp'))';
                
                if isa(mosaic,'rgcSubunit');
                    arbitraryScaleFactor = 73.73;
                else
                    arbitraryScaleFactor = 1;
                end
                
                % MAKE ZERO MEAN, maybe get rid of this
                fullResponseRSCombined = arbitraryScaleFactor*(fullResponseRSCombined - repmat(mean(fullResponseRSCombined,2),1,size(fullResponseRSCombined,2)));
               
                
%                 fullResponseRSCombinedRange = arbitraryScaleFactor*max(fullResponseRSCombined(:) - min(fullResponseRSCombined(:)));
%                 fullResponseRSCombined = (fullResponseRSCombinedRange/2)*fullResponseRSCombined;
                
                % figure; plot(1*fullResponseRSCombined(1328,:),'linewidth',3); hold on; plot(fullResponseRSCenter(1328,:),':r','linewidth',3)

                % Specify starting and ending time coordinates
                startPoint = 1; endPoint = nSamples;%+length(temporalIRCenter)-1;
                fullResponseRSRGB = zeros(spResponseSize(1)*spResponseSize(2),length(startPoint:endPoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCombined(:,startPoint:endPoint);
                                
            else            
                warning('Not circular conv, see fullConvolve.m');
                % if the temporal impulse responses for center and surround are different, convolve both               
% 
%                 K  = fittedGLM.linearfilters.Stimulus.Filter;
%                 
%                 K  = reshape(K, [ROI_pixels, length(frame_shifts)]);
                
%                 KX = zeros(13*13, length(spResponseCenterRS));
%                 for i_pixel = 1:13*13
%                     X_frame_shift = prep_timeshift(X_frame(i_pixel,:),0:29);
%                     tfilt = K(i_pixel,:);
% %                     KX(i_pixel,:) = tfilt * X_frame_shift;
%                     KX(i_pixel,:) = tfilt * temporalIRCenter;
%                 end
%                 lcif_kx_frame = sum(KX,1);
% % % % % %  Works with conv                
                fullResponseRSCenter = convn(spResponseCenterRS, temporalIRCenter','full');
                fullResponseRSSurround = convn(spResponseSurroundRS, temporalIRSurround','full');

%                 fullResponseRSCenter = ifft(fft(spResponseCenterRS).*repmat(fft(temporalIRCenter,169),1,1101));
%                 fullResponseRSSurround = ifft(fft(spResponseSurroundRS).*repmat(fft(temporalIRSurround,169),1,1101));
                
% % % % % Works with fft
%                 spResponseCenterRSp = [spResponseCenterRS zeros([size(spResponseCenterRS,1) size(temporalIRCenter,1)])];
%                 temporalIRCenterp = repmat([temporalIRCenter' zeros([size(temporalIRCenter,2) size(spResponseCenterRS,2)])],169,1);
%                 
%                 fullResponseRSCenter = ifft(fft(spResponseCenterRSp').*fft(temporalIRCenterp'))';
%                 
%                 spResponseSurroundRSp = [spResponseSurroundRS zeros([size(spResponseSurroundRS,1) size(temporalIRSurround,1)])];
%                 temporalIRSurroundp = repmat([temporalIRSurround' zeros([size(temporalIRSurround,2) size(spResponseSurroundRS,2)])],169,1);
%                 
%                 fullResponseRSSurround = ifft(fft(spResponseSurroundRSp).*fft(temporalIRSurroundp));
% % % % % % %                
                
                % Specify starting and ending time coordinates
%                 startPoint = length(temporalIRCenter)-1; endPoint = nSamples+length(temporalIRCenter)-1;
                startPoint = 1; endPoint = nSamples;%+length(temporalIRCenter)-1;
                % Take difference between center and surround response                
                fullResponseRSRGB = zeros(spResponseSize(1)*spResponseSize(2),length(startPoint:endPoint));
                fullResponseRSRGB(:,:,rgbIndex) = fullResponseRSCenter(:,startPoint:endPoint) - fullResponseRSSurround(:,startPoint:endPoint);
                
            end            
        end      
                        
        if isa(mosaic,'rgcPhys')
            
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
            % fullResponse{xcell,ycell,1} = mean(fullResponseRS);       % this is the only difference from the elseif block
            
            % Need to normalize by 13x13 RF size in lab physio code
            % fullResponse{xcell,ycell,1} = (13^2/size(fullResponseRS,1))*sum(fullResponseRS) + mosaic.tonicDrive{xcell,ycell};   
            fullResponse{xcell,ycell,1} = sum(fullResponseRS) + mosaic.tonicDrive{xcell,ycell};   
            % % fullResponse for RGB
            fullResponse{xcell,ycell,2} =  reshape(fullResponseRSRGB, spResponseSize(1), spResponseSize(2), size(fullResponseRSRGB,2), size(fullResponseRSRGB,3));
            if ~isa(mosaic, 'rgcLinear'); 
                genFunction = mosaicGet(mosaic, 'generatorFunction');
                nlResponse{xcell,ycell} = genFunction(mean(fullResponseRS,1));
            else
                
                nlResponse{xcell,ycell} = 0;
            end
        end
        
       
                                      
    end
end
toc