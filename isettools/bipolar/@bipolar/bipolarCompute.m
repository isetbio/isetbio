function obj = bipolarCompute(obj, os, varargin)
% Compute bipolar responses
% 
%   Still under active development
%   Should become:  bipolar.compute(coneMosaic,varargin);
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
p.addRequired('os', @(x) isa(x, 'outerSegment'));  % Should be coneMosaic

% parse
p.parse(obj, os, varargin{:});

%% Spatial filtering and subsampling
% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

% Get zero mean cone current signal
osSig = RGB2XWFormat(os.coneCurrentSignal);
osSigRSZM = bsxfun(@minus, osSig, mean(osSig, 2));
osSigZM = reshape(osSig, size(os.coneCurrentSignal));

%% Map cone positions to appropriate centers or surrounds of RGC RFs

% For offDiffuse, onDiffuse and onMidget, remove S cone inputs by replacing
% with nearest L/M input. For offMidget, keep S cones but scale down by
% 75%. For onSBC, only S cone inputs to center, only L/M cone inputs to
% surround.

switch obj.cellType
    % Remove S cone input for these types of bipolars
    case{'offDiffuse','onDiffuse','onMidget'}
        
        lmConeIndices = find(obj.coneMosaic ==2 | obj.coneMosaic == 3);
        sConeIndices = find(obj.coneMosaic==4);
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMSurround   = osSigRSZM;        
        
        [rLM,cLM]=ind2sub(size(obj.coneMosaic),lmConeIndices);
        [rS,cS]=ind2sub(size(obj.coneMosaic),sConeIndices);
        
        % Set center and surround to only have LM cones
        lmConeDist = sqrt((repmat(rLM,[1 length(rS)]) - repmat(rS',[length(rLM) 1])).^2 - (repmat(cLM,[1 length(cS)]) - repmat(cS',[length(cLM) 1])).^2);
                
        [mindlm,minindlm] = min(lmConeDist);
        % Plot S cones mapped to LM cones
        % vcNewGraphWin; scatter(rS(:),cS(:),20,lmConeIndices(minindlm),'filled')
        % colormap([rand(length(lmConeIndices),3)])
        osSigRSZMCenter(sConeIndices,:) = osSigRSZMCenter(lmConeIndices(minindlm),:);
        osSigRSZMSurround(sConeIndices,:) = osSigRSZMSurround(lmConeIndices(minindlm),:);
    % Keep S cone input for off Midget but only weight by 0.25
    case{'offMidget'}
        sConeIndices = find(obj.coneMosaic==4);
        minval = min(osSigRSZM(:));
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMCenter(sConeIndices,:)   = 0.25*(osSigRSZMCenter(sConeIndices,:)-minval)+minval;
        
        osSigRSZMSurround   = osSigRSZM;
        osSigRSZMSurround(sConeIndices,:) = 0.25*(osSigRSZMCenter(sConeIndices,:)-minval)+minval;
    % Make nearest S cones the center for SBCs, only L and M cones in
    % surround
    case{'onSBC'}        
        lmConeIndices = find(obj.coneMosaic ==2 | obj.coneMosaic == 3);
        sConeIndices = find(obj.coneMosaic==4);
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMSurround   = osSigRSZM;        
        
        % Set center to only have S cones
        
        [rLM,cLM]=ind2sub(size(obj.coneMosaic),lmConeIndices);
        [rS,cS]=ind2sub(size(obj.coneMosaic),sConeIndices);
        
        for sind = 1:length(sConeIndices)
            sConeDist(sind,:) = sqrt((rLM - rS(sind)).^2 + (cLM - cS(sind)).^2);
        end
    
        [mind,minind] = min(sConeDist);
        % Plot LM cones mapped to S cones
        % vcNewGraphWin; scatter(rLM(:),cLM(:),20,sConeIndices(minind),'filled')
        % colormap([rand(length(sConeIndices),3)])
        osSigRSZMCenter(lmConeIndices,:) = osSigRSZMCenter(sConeIndices(minind),:);
        
        % Set surround to only have LM cones
        lmConeDist = sqrt((repmat(rLM,[1 length(rS)]) - repmat(rS',[length(rLM) 1])).^2 - (repmat(cLM,[1 length(cS)]) - repmat(cS',[length(cLM) 1])).^2);
        
        [mindlm,minindlm] = min(lmConeDist);
        % Plot S cones mapped to LM cones
        % vcNewGraphWin; scatter(rS(:),cS(:),20,lmConeIndices(minindlm),'filled')
        % colormap([rand(length(lmConeIndices),3)])
        osSigRSZMSurround(sConeIndices,:) = osSigRSZMSurround(lmConeIndices(minindlm),:);
        
end

osSigZMCenter = reshape(osSigRSZMCenter,size(os.coneCurrentSignal));
osSigZMSurround = reshape(osSigRSZMSurround,size(os.coneCurrentSignal));

% Convolve spatially every frame
spatialResponseCenter = ieSpaceTimeFilter(osSigZMCenter, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(osSigZMSurround, obj.sRFsurround);

% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, strideSubsample);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, strideSubsample);

%% Temporal filtering

% Reshape for temporal convolution
szSubSample = size(spatialSubsampleCenter);

if numel(szSubSample)<3; szSubSample(3) = 1; end;
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

spatialSubsampleCenterRS = [repmat(spatialSubsampleCenterRS(:,1),1,1).*ones(size(spatialSubsampleCenterRS,1),1) spatialSubsampleCenterRS];
spatialSubsampleSurroundRS = [repmat(spatialSubsampleSurroundRS(:,1),1,1).*ones(size(spatialSubsampleSurroundRS,1),1) spatialSubsampleSurroundRS];    

switch obj.filterType
    case 1        
        % RDT initialization
        rdt = RdtClient('isetbio');
        rdt.crp('resources/data/rgc');
        if strcmpi(obj.cellType,'offDiffuse')
            data = load([isetRootPath '/data/bipolar/bipolarFilt_200_OFFP_2013_08_19_6_all.mat']);
        else
            data = load([isetRootPath '/data/bipolar/bipolarFilt_200_ONP_2013_08_19_6_all.mat']);
        end
        bipolarFiltMat = data.bipolarFiltMat;
        
        switch ieParamFormat(obj.cellType)
            case {'offdiffuse','offmidget'}
                bipolarFilt = -mean(bipolarFiltMat)';
            case {'ondiffuse','onmidget','onsbc'}
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
% bipolarOutputCenterRSLong = ifft(fft(bipolarOutputCenterRSLongZP').*fft(bipolarFiltZP'))';
% bipolarOutputSurroundRSLong = ifft(fft(bipolarOutputSurroundRSLongZP').*fft(bipolarFiltZP'))';
% 
% bipolarOutputCenterRS = bipolarOutputCenterRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);
% bipolarOutputSurroundRS = bipolarOutputSurroundRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);

% % % % % 
bipolarOutputCenterRS = convn(bipolarFilt',spatialSubsampleCenterRS,'same');
bipolarOutputSurroundRS = convn(bipolarFilt',spatialSubsampleSurroundRS,'same');

% bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','same');
% bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','same');
if size(spatialSubsampleCenterRS,2) < size(bipolarFilt,1)
    
    bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','full');
    bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','full');
    
    bipolarOutputCenterRS = bipolarOutputCenterRS(:,floor(size(bipolarFilt,1)/2):end);
    bipolarOutputSurroundRS = bipolarOutputSurroundRS(:,floor(size(bipolarFilt,1)/2):end);
    
elseif size(bipolarOutputCenterRS,2) > floor(size(bipolarFilt,1)/2)
    bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','same');
    bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','same');

    bipolarOutputCenterRS = bipolarOutputCenterRS(:,1:end-floor(size(bipolarFilt,1)/2));
    bipolarOutputSurroundRS = bipolarOutputSurroundRS(:,1:end-floor(size(bipolarFilt,1)/2));

end
% % % % % % 

% Rezero
bipolarOutputCenterRSRZ = ((bipolarOutputCenterRS-repmat(mean(bipolarOutputCenterRS,2),1,size(bipolarOutputCenterRS,2))));
bipolarOutputSurroundRSRZ = ((bipolarOutputSurroundRS-repmat(mean(bipolarOutputSurroundRS,2),1,size(bipolarOutputSurroundRS,2))));

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));

%% Calculate contrast gain adjustment

%% Attach output to object

obj.responseCenter = obj.rectificationCenter(bipolarOutputLinearCenter);
obj.responseSurround = obj.rectificationSurround(bipolarOutputLinearSurround);

end