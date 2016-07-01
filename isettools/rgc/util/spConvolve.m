function  [responseC, responseS] = spConvolve(rgcMosaic, stimulus, bipolarS)
% 2D spatial convolution for the receptive field
% 
%   [responseS, responseS] = spConvolve(rgcMosaic, stimulus, [bipolarSurround]);
%   
% Convolves that 2D stimulus with the spatial center and surround of each
% cell in the rgcMosaic.  This is applied independently to each of the
% temporal frames in the stimulus.
% 
% Inputs: 
%   rgcMosaic:  RGC ganglion cell mosaic (@rgcmosaic)
%   stimulus:   (x,y,t) usually the bipolar center current, but this can be
%               another stimulus (e.g., the outer segment current, in
%               which case there is no bipolarS
%   bipolarS:   (x,y,t) bipolar surround current (optional)
% 
% Outputs: 
%   responseC:  Mosaic of center responses over time
%   responseS:  Mosaic of surround responses over time
% 
% Example:
%   
%   [rC, rS] = spConvolve(rgc.mosaic{1}, cMosaic.os.coneCurrentSignal);
%   [rC, rS] = spConvolve(rgc.mosaic{1}, bp.responseCenter, bp.responseSurround);
%     
% (c) isetbio
% 09/2015 JRG

nSamples = size(stimulus,3);
channelSize = size(stimulus,4); % if RGB, channelSize = 3

nCells = size(rgcMosaic.cellLocation);

fprintf('     \n');
fprintf('Spatial Convolution, %s:     \n', rgcMosaic.cellType);

responseC = cell(nCells(1),nCells(2));
responseS = cell(nCells(1),nCells(2));
for rgbIndex = 1:channelSize 
    tic
    
    % Calculating center and response at each cell
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            
            % Initialize cell arrays
            responseC{xcell,ycell} = zeros([size(rgcMosaic.sRFcenter{xcell,ycell}) nSamples channelSize]); 
            responseS{xcell,ycell} = zeros([size(rgcMosaic.sRFcenter{xcell,ycell}) nSamples channelSize]); 
            
            % Loop over each time step
            for tt = 1:nSamples
                
                % Get RF weights for this cell
                spRFcenter   = rgcMosaic.sRFcenter{xcell,ycell};
                spRFsurround = rgcMosaic.sRFsurround{xcell,ycell};
                
                % Get cell locations in stimulus referred coordinate frame
                [stimX, stimY, offset] = rgcMosaic.stimPositions(xcell,ycell);
                
                % Check rgc inputs are within size of stimulus
                gz = find((stimX-offset(1))>=1 & ...
                          (stimY-offset(2))>=1 & ...
                          (stimX-offset(1))<=size(stimulus,1) & ...
                          (stimY-offset(2))<=size(stimulus,2) );
                
                % Extract 2D image
                spStim = squeeze(stimulus(floor(stimX(gz)-offset(1)),floor(stimY(gz)-offset(2)),tt,rgbIndex));
                
                % Extract 2D image for surround computation, same as center
                if isa(rgcMosaic, 'rgcSubunit') && exist('spStimSurr','var')
                    spStimSurr = squeeze(bipolarS(floor(stimX(gz)-offset(1)),floor(stimY(gz)-offset(2)),tt,rgbIndex));
                end
                
                % Do elementwise product of spatial RF with stimulus frame
                % of the same size
               
                if isa(rgcMosaic, 'rgcPhys')
                    % Data driven case
                    responseC{xcell,ycell}(gz,gz,tt,rgbIndex) = spRFcenter(gz,gz).*spStim;
                    
                elseif isa(rgcMosaic, 'rgcSubunit')
                    % For computation with a subunit rgc, apply recitfying nonlinearities
                
                    spRC = (spRFcenter(gz,gz).*(spStim-1*mean(spStim(:))));
                    
                    if exist('spStimSurr','var')
                        spRS = (spRFsurround(gz,gz).*-0*(spStimSurr-1*mean(spStimSurr(:))));
                        responseC{xcell,ycell}(gz,gz,tt,rgbIndex) = rgcMosaic.rectifyFunction(sum(spRC(:)))./length(gz)^2;
                        responseS{xcell,ycell}(gz,gz,tt,rgbIndex) = rgcMosaic.rectifyFunction(sum(spRS(:)))./length(gz)^2;

                    else 
                        spRS = (spRFsurround(gz,gz).*-(spStim-1*mean(spStim(:))));
                        responseC{xcell,ycell}(1,1,tt,rgbIndex) = 10*sum(rgcMosaic.rectifyFunction(spRC(:)))./length(gz)^2;
                        responseS{xcell,ycell}(1,1,tt,rgbIndex) = 10.*sum(rgcMosaic.rectifyFunction(spRS(:)))./length(gz)^2;
                    end
                else
                    % Typical linear, LNP and GLM cases
                    responseC{xcell,ycell}(:,:,tt,rgbIndex) = conv2(spRFcenter, spStim-1*mean(spStim(:)), 'same');
                    responseS{xcell,ycell}(:,:,tt,rgbIndex) = conv2(spRFsurround, spStim-1*mean(spStim(:)), 'same');

                end
            end
        end
    end
    toc
end
