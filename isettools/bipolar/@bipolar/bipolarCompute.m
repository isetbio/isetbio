function obj = bipolarCompute(obj, os, varargin)
% Computes the responses of the bipolar object. 
% 
% The outersegment input contains frames of cone mosaic signal at a
% particular time step. The bipolar response is found by first convolving
% the center and surround Gaussian spatial receptive fields of the bipolar
% cell within each cone signal frame. Then, the resulting signal is put
% through the weighted temporal differentiator in order to result in an
% impulse response that approximates the IR of the RGC.
% 
% Particular options that could be employed are rezeroing of the signal at
% the end of the temporal computation as well as rectification on the
% output signal.
% 
% 5/2016 JRG (c) isetbio team

%% parse input parameters
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'bipolar'));
p.addRequired('os', @(x) isa(x, 'outerSegment'));

% parse
p.parse(obj, os, varargin{:});

%% Spatial filtering and subsampling
% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

% Get zero mean cone current signal
osSig = RGB2XWFormat(os.coneCurrentSignal);
osSigRSZM = bsxfun(@minus, osSig, mean(osSig, 2));
osSigZM = reshape(osSig, size(os.coneCurrentSignal));

%% Spatial response

%%%%%% Old can be removed
% % Get zero mean cone current signal
% osSigRS = reshape(os.coneCurrentSignal, size(os.coneCurrentSignal,1)*size(os.coneCurrentSignal,2),size(os.coneCurrentSignal,3));
% osSigRSZM = osSigRS - repmat(mean(osSigRS,2),1,size(osSigRS,2));
%%%%%%%%%

% osSigZM = reshape(osSigRSZM,size(os.coneCurrentSignal));


% Take product of outer segment current with the cone mask in order to
% scale contributions of S cones appropriately. For off midget cells, scale
% cone current by a factor of 0.25. For SBCs
% scaledCurrent = bp.coneMask

% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

switch obj.cellType
    case{'offDiffuse','onDiffuse','onMidget'}
        lmConeIndices = find(obj.coneMosaic ==2 | obj.coneMosaic == 3);
        sConeIndices = find(obj.coneMosaic==4);

        osSigRSZMCenter   = osSigRSZM;
        % osSigRSZMCenter(sConeIndices) = mean(osSigRSZMCenter(lmConeIndices));
        
        osSigRSZMSurround   = osSigRSZM;
        % osSigRSZMSurround(sConeIndices) = mean(osSigRSZMSurround(lmConeIndices)        
             
        currentSignalCenter = osSigRSZMCenter;
        currentSignalSurround = osSigRSZMSurround;
        
        sRepeat = find(diff(sConeIndices)==1);
        
        for fr = 1:size(currentSignalCenter,2)
            coneCurrentFrame = osSigRSZMCenter(:,fr);
            % meanLM = mean(coneCurrentFrame(lmConeIndices));
            coneCurrentCenter = coneCurrentFrame;
            coneCurrentSurround = coneCurrentFrame;
            coneCurrentCenter(sConeIndices) = coneCurrentFrame(sConeIndices-1);
            coneCurrentCenter(sConeIndices(sRepeat)+1) = coneCurrentFrame(sConeIndices(sRepeat+1)-2);
            coneCurrentSurround(sConeIndices) = coneCurrentFrame(sConeIndices-1);
            coneCurrentSurround(sConeIndices(sRepeat+1)) = coneCurrentFrame(sConeIndices(sRepeat+1)-2);
            osSigRSZMCenter(:,fr) = coneCurrentCenter;
            osSigRSZMSurround(:,fr) = coneCurrentSurround;
        end

        
    case{'offMidget'}
        sConeIndices = find(obj.coneMosaic==4);
        minval = min(osSigRSZM(:));
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMCenter(sConeIndices,:)   = 0.25*(osSigRSZMCenter(sConeIndices,:)-minval)+minval;
        
        osSigRSZMSurround   = osSigRSZM;
        osSigRSZMSurround(sConeIndices,:) = 0.25*(osSigRSZMCenter(sConeIndices,:)-minval)+minval;
    
    case{'onDiffuseSBC'}        
        lmConeIndices = find(obj.coneMosaic ==2 | obj.coneMosaic == 3);
        sConeIndices = find(obj.coneMosaic==4);
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMSurround   = osSigRSZM;
        
        currentSignalCenter = osSigRSZMCenter;
        currentSignalSurround = osSigRSZMSurround;
        
        sRepeat = find(diff(sConeIndices)==1);
        
        for fr = 1:size(currentSignalCenter,2)
            coneCurrentFrame = osSigRSZMCenter(:,fr);
            % meanLM = mean(coneCurrentFrame(lmConeIndices));
%             coneCurrentCenter = coneCurrentFrame;
            coneCurrentSurround = coneCurrentFrame;
%             coneCurrentCenter(sConeIndices) = coneCurrentFrame(sConeIndices-1);
%             coneCurrentCenter(sConeIndices(sRepeat)+1) = coneCurrentFrame(sConeIndices(sRepeat+1)-2);
            coneCurrentSurround(sConeIndices) = coneCurrentFrame(sConeIndices-1);
            coneCurrentSurround(sConeIndices(sRepeat+1)) = coneCurrentFrame(sConeIndices(sRepeat+1)-2);
%             osSigRSZMCenter(:,fr) = coneCurrentCenter;
            osSigRSZMSurround(:,fr) = coneCurrentSurround;
        end
        
        [rLM,cLM]=ind2sub(size(obj.coneMosaic),lmConeIndices);
        [rS,cS]=ind2sub(size(obj.coneMosaic),sConeIndices);
        
        for sind = 1:length(sConeIndices)
            sConeDist(sind,:) = sqrt((rLM - rS(sind)).^2 + (cLM - cS(sind)).^2);
        end
    
        [mind,minind] = min(sConeDist);
        % vcNewGraphWin; scatter(rLM(:),cLM(:),20,sConeIndices(minind),'filled')
        % colormap([rand(length(sConeIndices),3)])
        osSigRSZMCenter(lmConeIndices,:) = osSigRSZMCenter(sConeIndices(minind),:);
end

osSigZMCenter = reshape(osSigRSZMCenter,size(os.coneCurrentSignal));
osSigZMSurround = reshape(osSigRSZMSurround,size(os.coneCurrentSignal));

% osSigZM = reshape(osSigRSZM,size(os.coneCurrentSignal));
% osSigZMCenter = osSigZM;
% osSigZMSurround = osSigZM;

spatialResponseCenter = ieSpaceTimeFilter(osSigZMCenter, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(osSigZMSurround, obj.sRFsurround);

% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, strideSubsample);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, strideSubsample);

%% Temporal filtering
% Apply the weighted differentiator to the output of the spatial response

% Reshape for temporal convolution

%%% Implement if works
% <<<<<<< HEAD
% % Compute spatial center/surround response
% sCenter = convn(osSigZM, obj.sRFcenter, 'same');
% sSurround = convn(osSigZM, obj.sRFsurround, 'same');
% 
% % Subsample to pull out individual bipolars
% stride = size(obj.sRFcenter, 1);
% sCenter = sCenter(1:stride:size(sCenter, 1), 1:stride:size(sCenter, 2), :);
% sSurround = sSurround(1:stride:size(sSurround, 1), ...
%     1:stride:size(sSurround, 2), :);

% [sCenter, r, c] = RGB2XWFormat(sCenter);
% sSurround = RGB2XWFormat(sSurround);
% 
% % Repeat input matrix for filtering
% sCenter = [repmat(sCenter(:,1), 1, 1) sCenter];
% sSurround = [repmat(sSurround(:,1), 1, 1) sSurround];    
% 
% % load filters
% if obj.filterType == 1  % average filter from measurement data
%     if strcmpi(obj.cellType, 'offDiffuse')
%         data = load('bipolarFilt.mat', 'bipolarOFFP');
%         bipolarFilt = -data.bipolarOFFP(:)';
%     else
%         data = load('bipolarFilt.mat', 'bipolarONP');
%         bipolarFilt = -data.bipolarONP(:)';
%     end
%     
%     % Interpolate filter from data to get correct sample rate
%     originalTimeStep = 0.008;
%     bpLength = length(bipolarFilt);
%     bipolarFilt = interp1(originalTimeStep:originalTimeStep:originalTimeStep*bpLength,bipolarFilt,obj.timeStep:obj.timeStep:originalTimeStep*bpLength);
%     bipolarFilt(isnan(bipolarFilt)) = 0;
% =======


szSubSample = size(spatialSubsampleCenter);
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

% obj.temporalDelay = 0;
temporalDelay = 0;
% Zero pad to allow for delay
spatialSubsampleCenterRS = [repmat(spatialSubsampleCenterRS(:,1),1,(1e-3/os.timeStep)*temporalDelay + 1).*ones(size(spatialSubsampleCenterRS,1),(1e-3/os.timeStep)*temporalDelay + 1) spatialSubsampleCenterRS];
spatialSubsampleSurroundRS = [repmat(spatialSubsampleSurroundRS(:,1),1,(1e-3/os.timeStep)*temporalDelay + 1).*ones(size(spatialSubsampleSurroundRS,1),(1e-3/os.timeStep)*temporalDelay + 1) spatialSubsampleSurroundRS];    

switch obj.filterType
    case 1
        
        % RDT initialization
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        if strcmpi(obj.cellType,'offDiffuse')
            data = rdt.readArtifact('bipolarFilt_200_OFFP_2013_08_19_6_all', 'type', 'mat');            
            
        else
            data = rdt.readArtifact('bipolarFilt_200_ONP_2013_08_19_6_all', 'type', 'mat');            
            
        end
        bipolarFiltMat = data.bipolarFiltMat;
        
        switch ieParamFormat(obj.cellType)
            case {'offdiffuse','offmidget'}
                bipolarFilt = -mean(bipolarFiltMat)';
            case {'ondiffuse','onmidget','ondifusesbc'}
                bipolarFilt = -mean(bipolarFiltMat)';
            otherwise
                error('Unknown bipolar cell type');
        end
        
        
    case  2
        load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/irGLM.mat');
        if strcmpi(obj.cellType, 'offDiffuse')
            bipolarFilt = irGLM;
        else
            bipolarFilt = -irGLM;
        end

    case 3
        % RDT initialization
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        if strcmpi(obj.cellType,'offDiffuse')
            data = rdt.readArtifact('bipolarFilt_200_OFFP_2013_08_19_6_all', 'type', 'mat');
        else
            data = rdt.readArtifact('bipolarFilt_200_ONP_2013_08_19_6_all', 'type', 'mat');
        end
        bipolarFiltMat = data.bipolarFiltMat;
        load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/bipolarFilt_200_OFFP_2013_08_19_6_all_linear.mat');
   
    
    case 4  % sampled at 150 fr/sec for impulse response
    
        data = load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/bipolarFilt_200_ONP_2013_08_19_6_all_linear_fr150.mat');
        bipolarFilt = (data.bipolarFiltMat(obj.cellLocation,:)');
end


% bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','same');
% bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','same');
% bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','full');
% bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','full');

% % temporal filtering
% tCenter = conv2(sCenter, bipolarFilt, 'same');
% tSurround = conv2(sSurround, bipolarFilt, 'same');
% 
% % zero centering the output
% tCenter = reshape(bsxfun(@minus, tCenter, mean(tCenter, 2)), r, c, []);
% tSurround = reshape(bsxfun(@minus,tSurround,mean(tSurround, 2)), r, c, []);

% %% Rectification
% obj.responseCenter = obj.rectificationCenter(tCenter);
% obj.responseSurround = obj.rectificationSurround(tSurround);
% 
% end

% bipolarFilt = (bipolarFiltMat(1,:)');
if size(spatialSubsampleCenterRS,2) > size(bipolarFilt,1)
    bipolarOutputCenterRSLongZP = [spatialSubsampleCenterRS];% zeros([size(spatialSubsampleCenterRS,1) size(bipolarFilt,1)])];
    bipolarOutputSurroundRSLongZP = [spatialSubsampleSurroundRS];% zeros([size(spatialSubsampleSurroundRS,1)-size(bipolarFilt,1)])];
    bipolarFiltZP = repmat([bipolarFilt; zeros([-size(bipolarFilt,1)+size(spatialSubsampleCenterRS,2)],1)]',size(spatialSubsampleCenterRS,1) ,1);
else

    bipolarOutputCenterRSLongZP = ([spatialSubsampleCenterRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleCenterRS,2)],1)',size(spatialSubsampleCenterRS,1),1)]);
    
    bipolarOutputSurroundRSLongZP = ([spatialSubsampleSurroundRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleSurroundRS,2)],1)',size(spatialSubsampleSurroundRS,1),1)]);
    bipolarFiltZP = repmat(bipolarFilt',size(spatialSubsampleSurroundRS,1),1);
    
end

% 
bipolarOutputCenterRSLong = ifft(fft(bipolarOutputCenterRSLongZP').*fft(bipolarFiltZP'))';
bipolarOutputSurroundRSLong = ifft(fft(bipolarOutputSurroundRSLongZP').*fft(bipolarFiltZP'))';

bipolarOutputCenterRS = bipolarOutputCenterRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);
bipolarOutputSurroundRS = bipolarOutputSurroundRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);


% Rezero
bipolarOutputCenterRSRZ = ((bipolarOutputCenterRS-repmat(mean(bipolarOutputCenterRS,2),1,size(bipolarOutputCenterRS,2))));
bipolarOutputSurroundRSRZ = ((bipolarOutputSurroundRS-repmat(mean(bipolarOutputSurroundRS,2),1,size(bipolarOutputSurroundRS,2))));
% figure; plot(conv(spatialSubsampleRS(50,:),obj.tIR'))

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));
% figure; plot(squeeze(bipolarOutputLinear(20,20,:)));

%% Calculate contrast gain adjustment

%% Attach output to object
% obj.responseCenter = os.coneCurrentSignal;
% obj.responseSurround = zeros(size(os.coneCurrentSignal));

obj.responseCenter = obj.rectificationCenter(bipolarOutputLinearCenter);
obj.responseSurround = obj.rectificationSurround(bipolarOutputLinearSurround);

% % Bipolar rectification 
% obj.responseCenter = os.coneCurrentSignal;
% obj.responseSurround = zeros(size(os.coneCurrentSignal));
% 
% obj.responseCenter = (bipolarOutputLinearCenter);
% obj.responseSurround = zeros(size(bipolarOutputLinearSurround));

% bipolarOutputRectifiedCenter = bipolarOutputLinearCenter.*(bipolarOutputLinearCenter>0);
% bipolarOutputRectifiedSurround = zeros(size(-bipolarOutputLinearSurround.*(bipolarOutputLinearSurround<0)));
% 
% % bipolarOutputRectifiedCenter = 1*abs(bipolarOutputLinearCenter);
% % bipolarOutputRectifiedSurround = zeros(size(bipolarOutputLinearSurround));
% 
% obj.responseCenter = bipolarOutputRectifiedCenter;
% obj.responseSurround = bipolarOutputRectifiedSurround;

