function obj = bipolarCompute(obj, cmosaic, varargin)
% Compute bipolar responses
% 
%   Still under active development
%   Should become:  bipolar.compute(coneMosaic,varargin);
%
% Inputs:
%   obj:       a bipolar object
%   cmosaic:   coneMosaic  (N.B.  We allow an os for backward
%              compatibility, but that will be deprecated).
% 
% Anatomical connections:
%  The bipolar cells are classified into several types
%    
%  * on/off diffuse, which connect to parasol RGCs
%  * on/off midget, which connect to midget RGCs
%  * on small bistratified, which connect to S-cone bistratified
%
% This function forces the receptive field properties in terms of cone
% connections to match up correctly with the cone mosaic
%
% Computations:
%  The outersegment input contains frames of cone mosaic signal at a
%  particular time step. The bipolar response is found by first convolving
%  the center and surround Gaussian spatial receptive fields of the bipolar
%  cell within each cone signal frame. Then, the resulting signal is put
%  through the weighted temporal differentiator in order to result in an
%  impulse response that approximates the IR of the RGC.
% 
% TODO:
% Particular options that could be employed are rezeroing of the signal at
% the end of the temporal computation as well as rectification on the
% output signal.
% 
% 5/2016 JRG (c) isetbio team

%% parse input parameters
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'bipolar'));
p.addRequired('cmosaic', @(x) isa(x, 'outerSegment') | isa(x, 'coneMosaic'));  

% parse
p.parse(obj, cmosaic, varargin{:});

% The input object should be coneMosaic, but it can also be an OS for
% backwards compatibility for now.
if isa(cmosaic,'coneMosaic'),     os = cmosaic.os;
else                              os = cmosaic;
end

%% Spatial filtering and subsampling
% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

% Zero-mean the cone current signal at each cone

% This places the cone 3D matrix into a coneNumber x time matrix
osSig = RGB2XWFormat(os.coneCurrentSignal);

% Typically there 
if size(osSig,2) > 1
    % Typical case.  Substract the mean over time of each cone signal from
    % itself. 
    osSigRSZM = bsxfun(@minus, osSig, mean(osSig, 2));
else
    % Sometimes there is only one time point, so don't subtract it from
    % itself (which would be zero).
    osSigRSZM = osSig;
end

%% Enfoce anatomical requirements on cone connections

% Rules:
%
%  off Diffuse, on Diffuse and on Midget - remove S cone inputs
%    These are replaced with nearest L/M input. 
%  For offMidget, keep S cones but scale the connection strength down by 75%. 
%  For onSBC, only S cone inputs to center, only L/M cone inputs to surround.
%
% Citation:  See bipolar.m

% We need to let people change this.  It should be a choice when people run
% the code rather than buried in here.  It may matter a lot, or not, but it
% should be open to people experimenting with it.
switch obj.cellType
    case{'offDiffuse','onDiffuse','onMidget'}

        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMSurround = osSigRSZM;

        % Remove S cone input for these types of bipolars

        % Find the locations (row, col) of the different cone types
        [~,~,S] = coneTypeLocations(cmosaic,'val','index');
        
        % Zero the photocurrent of the S cones. Do this for both the center
        % and the surround.
        z = zeros(length(S),size(osSigRSZM,2));
        osSigRSZMCenter(S(:),:)   = z;
        osSigRSZMSurround(S(:),:) = z;
        
    case{'offMidget'}
        % Keep S cone input for off Midget but only weight by 0.25
        
        % Find the locations (row, col) of the different cone types
        [~,~,S] = coneTypeLocations(cmosaic,'val','index');
        
        minval = min(osSigRSZM(:));
        
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMCenter(S,:)   = 0.25*(osSigRSZMCenter(S,:)-minval)+minval;
        
        osSigRSZMSurround   = osSigRSZM;
        osSigRSZMSurround(S,:) = 0.25*(osSigRSZMSurround(S,:)-minval)+minval;

    case{'onSBC'}  
        % Set L and M cones to zero in SBC center, set S cones to zero in
        % SBC surround.
        % Find the locations (row, col) of the different cone types
        [L,M,S] = coneTypeLocations(cmosaic,'val','index');
        LM = [L; M];
        
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMSurround = osSigRSZM;        
        
        minval = min(osSigRSZM(:));
        % Set center to only have S cones
        
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMCenter(LM,:)   = minval*ones(size(osSigRSZMCenter(LM,:)));
        
        osSigRSZMSurround   = osSigRSZM;
        osSigRSZMSurround(S,:)   = minval*ones(size(osSigRSZMSurround(S,:)));
                       
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

%% Pull bipolar temporal filter from RDT
% The bipolar temporal filter is a result of the deconvolution of the
% linear cone temporal response and the linear RGC temporal response. There
% are several different bipolar filters that the user can select based on
% the type of simulation (simulated parameters or physiological
% parameters).
% 
% For the bipolar filter from simulated RGC parameters, the ideal RGC
% temporal response was generated using the code by Jonathan Pillow. For
% the bipolar filter from physiological RGC parameters, a bipolar filter
% for each individual RGC was generated. These may be averaged together or
% kept separate for certain computations (impulse response, comparison of
% RGC responses from isetbio to RGC responses from the Chichilnisky Lab's
% code.
% 
% In order to set the filterType property for the bipolar, the user must
% pass it as a parameter when the bipolar mosaic is created.

switch obj.filterType
    case 1
        % The basic physiology response case. Bipolar filters were
        % deconvolved from the measured temporal impulse response of each
        % cell in the mosaic and the linear cone temporal response. The
        % mean of the bipolar temporal filters for the whole mosaic is used
        % as the ideal bipolar filter.
        
        % There are different bipolar filters for on p/m and off p/m cells
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
        % The basic simulated response case. Bipolar filters were
        % deconvolved from the ideal temporal impulse response from the
        % Pillow simulation code.
        load([isetRootPath  '/data/bipolar/irGLM.mat']);
        if strcmpi(obj.cellType, 'offDiffuse')
            bipolarFilt = irGLM;
        else
            bipolarFilt = -irGLM;
        end

    case 3
        % The physiology response case that allows the individual bipolar
        % filters to be used. This case is for testing the impulse response
        % of the isetbio RGC and comparing to the impulse response measured
        % in physiology, or for matching RGC responses with code from the
        % Chichilnisky Lab.
        if strcmpi(obj.cellType,'offDiffuse')            
            data = load([isetRootPath '/data/bipolar/bipolarFilt_200_OFFP_2013_08_19_6_all.mat']);
        else
            data = load([isetRootPath '/data/bipolar/bipolarFilt_200_ONP_2013_08_19_6_all.mat']);
        end
        
        bipolarFilt = -(data.bipolarFiltMat(obj.cellLocation,:)');

    case 4  
        % The same as case 3, but at a higher sampling rate,sampled at 
        % 150 fr/sec for the impulse response calculation.
        data = load([isetRootPath  '/data/bipolar/bipolarFilt_200_ONP_2013_08_19_6_all_linear_fr150.mat']);
        bipolarFilt = (data.bipolarFiltMat(obj.cellLocation,:)');
end

% bipolarFilt = (bipolarFiltMat(1,:)');

%% Zero pad filter or signal
% Zero paddding to allow for computations with FFT and IFFT (circular convolution). 
% This is primarily for handling the case when the stimulus is only one frame.

if size(spatialSubsampleCenterRS,2) > size(bipolarFilt,1)
    % Stimulus input longer than bipolar filter   
    bipolarOutputCenterRSLongZP = [spatialSubsampleCenterRS];% zeros([size(spatialSubsampleCenterRS,1) size(bipolarFilt,1)])];
    bipolarOutputSurroundRSLongZP = [spatialSubsampleSurroundRS];% zeros([size(spatialSubsampleSurroundRS,1)-size(bipolarFilt,1)])];
    
    % Zero pad the bipolar filter and repmat in order to convolve with each cone
    bipolarFiltZP = repmat([bipolarFilt; zeros([-size(bipolarFilt,1)+size(spatialSubsampleCenterRS,2)],1)]',size(spatialSubsampleCenterRS,1) ,1);
else
    % Stimulus shorter than bipolar filter (likely one frmae)
    % Zero pad the stimulus
    bipolarOutputCenterRSLongZP = ([spatialSubsampleCenterRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleCenterRS,2)],1)',size(spatialSubsampleCenterRS,1),1)]);
    
    % Repmat the bipolar filter for each cone input
    bipolarOutputSurroundRSLongZP = ([spatialSubsampleSurroundRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleSurroundRS,2)],1)',size(spatialSubsampleSurroundRS,1),1)]);
    bipolarFiltZP = repmat(bipolarFilt',size(spatialSubsampleSurroundRS,1),1);    
end

%% Compute the temporal response of the bipolar mosaic

% % Compute with FFT (circular convolution, no transient responses). This
% % works for reproducing the responses from the Chichilnisky Lab code.
% bipolarOutputCenterRSLong = ifft(fft(bipolarOutputCenterRSLongZP').*fft(bipolarFiltZP'))';
% bipolarOutputSurroundRSLong = ifft(fft(bipolarOutputSurroundRSLongZP').*fft(bipolarFiltZP'))';
% bipolarOutputCenterRS = bipolarOutputCenterRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);
% bipolarOutputSurroundRS = bipolarOutputSurroundRSLong;%(:,1:end-(1e-3/os.timeStep)*temporalDelay);

% Compute with convn (includes transient response). This is important for
% handling the one-frame stimulus case. 
% % % % % 
bipolarOutputCenterRS = convn(bipolarFilt',spatialSubsampleCenterRS,'same');
bipolarOutputSurroundRS = convn(bipolarFilt',spatialSubsampleSurroundRS,'same');

% bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','same');
% bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','same');
if size(spatialSubsampleCenterRS,2) < size(bipolarFilt,1)
    % Stimulus input longer than bipolar filter temporal length
    bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','full');
    bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','full');
    
    % Get rid of the transient onset response
    bipolarOutputCenterRS = bipolarOutputCenterRS(:,floor(size(bipolarFilt,1)/2):end);
    bipolarOutputSurroundRS = bipolarOutputSurroundRS(:,floor(size(bipolarFilt,1)/2):end);
    
elseif size(bipolarOutputCenterRS,2) > floor(size(bipolarFilt,1)/2)
    % Bipolar filter temporal length longer than stimulus input, probably
    % one frame
    bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,bipolarFilt','same');
    bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,bipolarFilt','same');
    % Get rid of the transient onset response
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