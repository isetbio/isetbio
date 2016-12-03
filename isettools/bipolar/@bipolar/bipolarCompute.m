function obj = bipolarCompute(obj, cmosaic, varargin)
% Compute bipolar continuous current responses
% 
%    bipolar.compute(coneMosaic,varargin);
%
% The bipolars act as a spatial-temporal function that converts the cone
% photocurrent into bipolar current that is delivered to the retinal
% ganglion cells.
%
% Inputs:
%   obj:       a bipolar object
%   cmosaic:   coneMosaic 
% 
% Key parameters
%
%  Cell type - The bipolar cells are classified into several types
%    
%  * on/off diffuse, which connect to parasol RGCs
%  * on/off midget, which connect to midget RGCs
%  * on small bistratified, which connect to S-cone bistratified
%
%  Processing - linear or rectified
%
%  The cone mosaic current contains a time series of the photocurrent
%  (coneMosaic.current). The bipolar response performs a spatio-temporal
%  separable convolution (first spatial filtering, then temporal filtering)
%  of the cone current to create the bipolar current.
%
%  The principal decision is whether the bipolar transformation is linear
%  or includes a rectification.  This is controlled by the obj.rectifyType
%  parameter. The rectification happens at the end of the function, after 
%  both spatial and temporal filtering.
%  REFERENCE: Meister option
%  
%  The spatial filtering is followed by a temporal filter that is selected
%  in order to match the impulse response that expect to find at the RGC
%  input.  REFERENCE: Chichilnisky option
% 
% 5/2016 JRG (c) isetbio team

%% parse input parameters

p = inputParser;
p.addRequired('obj', @(x) (isa(x, 'bipolar')));
p.addRequired('cmosaic', @(x) (isa(x, 'coneMosaic')));  

% parse - no options at this opint
p.parse(obj, cmosaic, varargin{:});

%% Spatial filtering and subsampling

% Convolve spatial RFs across the photo current of the cones in the mosaic

% This places the cone 3D matrix into a coneNumber x time matrix
osSig = RGB2XWFormat(cmosaic.current);

%% Zero-mean the input signal.
%
% BW thinks that the receptive field should govern how we map the input
% current to the output current.   If the RF has a zero mean, then spatial
% mean -> 0. If the RF has a unit mean then spatial mean -> mean
%
if size(osSig,2) > 1
    % Typical case.  Substract the mean over time of each cone signal from
    % itself. 
    osSigRSZM = bsxfun(@minus, osSig, mean(osSig, 2));
end

%% Enfoce anatomical rules on cone to bipolar connections

% Anatomical rules:
%
%  off Diffuse, on Diffuse and on Midget - These receive no S cone input
%  offMidget - keep S cones but scale the connection strength down by 75% 
%  onSBC     - S cone inputs to center, only L/M cone inputs to surround
%
% Citations:  See bipolar.m.  Wiki page <>

% TODO:  In the future we should set this up as a structure that we use to
% implement the anatomical rules.  Let's send in a struct that defines the
% anatomical rules (e.g., aRules) with slots that implement the kind of
% stuff listed above.
%
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
                
        minval = min(osSigRSZM(:));
        % Set center to only have S cones
        
        osSigRSZMCenter   = osSigRSZM;
        osSigRSZMCenter(LM,:)   = minval*ones(size(osSigRSZMCenter(LM,:)));
        
        osSigRSZMSurround   = osSigRSZM;
        osSigRSZMSurround(S,:)   = minval*ones(size(osSigRSZMSurround(S,:)));
                       
end

% Put the data back into RGB format, like RGB2XW()
osSigZMCenter = reshape(osSigRSZMCenter,size(cmosaic.current));
osSigZMSurround = reshape(osSigRSZMSurround,size(cmosaic.current));

%% Spatial convolution

% Full spatial convolution for every frame
spatialResponseCenter   = ieSpaceTimeFilter(osSigZMCenter, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(osSigZMSurround, obj.sRFsurround);

% Subsample in space to the resolution we expect for this bipolar mosaic
% The spacing is equal to the number of pixels that make up the center of
% the spatial receptive field.  This could be a settable parameter for
% others to experiment with, too.
spacing = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, spacing);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, spacing);

%% Temporal filtering

% Reshape for temporal convolution
szSubSample = size(spatialSubsampleCenter);

if numel(szSubSample)<3; szSubSample(3) = 1; end;
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

spatialSubsampleCenterRS = [repmat(spatialSubsampleCenterRS(:,1),1,1).*ones(size(spatialSubsampleCenterRS,1),1) spatialSubsampleCenterRS];
spatialSubsampleSurroundRS = [repmat(spatialSubsampleSurroundRS(:,1),1,1).*ones(size(spatialSubsampleSurroundRS,1),1) spatialSubsampleSurroundRS];    

%% Load bipolar temporal filter

% TODO:  Do the temporal interpolation of the bipolar time filter

% Bipolar filters were deconvolved from the measured temporal impulse response of each
% cell in the mosaic and the linear cone temporal response. The
% mean of the bipolar temporal filters for the whole mosaic is used
% as the ideal bipolar filter.
% 
% There are different bipolar filters for on p/m and off p/m cells
if strcmpi(obj.cellType,'offDiffuse')
    % Off parasol (off diffuse) only
    data = load([isetRootPath '/data/bipolar/bipolarFilt_200_OFFP_2013_08_19_6_all.mat']);
else
    % On parasol and all others
    data = load([isetRootPath '/data/bipolar/bipolarFilt_200_ONP_2013_08_19_6_all.mat']);
end

% The bipolarFiltMat is 100 different filters.  Fix that in the
% representation.  The script that builds this is about to be uploaded to
% the repository by JRG.
% bipolarFilt = -mean(data.bipolarFiltMat)';
% TODO: Remove bipolar.filterType parameter from object

bipolarFilt = bipolarFilter(obj, cmosaic);
warning('filter created here, see bipolarFilter.m; looks noncausal?');
%% Compute the temporal response of the bipolar mosaic

% Compute with convn (includes transient response). This is important for
% handling the one-frame stimulus case. 

% % Could problem be ordering of filter and input signal? Doesn't seem like
% % it.
% bipolarOutputCenterRSlong   = convn(spatialSubsampleCenterRS,  bipolarFilt','full');
% bipolarOutputSurroundRSlong = convn(spatialSubsampleSurroundRS,bipolarFilt','full');
    
bipolarOutputCenterRSlong   = convn(bipolarFilt',spatialSubsampleCenterRS,  'full');
bipolarOutputSurroundRSlong = convn(bipolarFilt',spatialSubsampleSurroundRS,'full');

warning('some sort of wrapping around end of time axis...');
bipolarOutputCenterRS   = bipolarOutputCenterRSlong(:,1:size(spatialSubsampleCenterRS,2));
bipolarOutputSurroundRS = bipolarOutputSurroundRSlong(:,1:size(spatialSubsampleCenterRS,2));
%% Format data

% Rezero
bipolarOutputCenterRSRZ = ((bipolarOutputCenterRS-repmat(mean(bipolarOutputCenterRS,2),1,size(bipolarOutputCenterRS,2))));
bipolarOutputSurroundRSRZ = ((bipolarOutputSurroundRS-repmat(mean(bipolarOutputSurroundRS,2),1,size(bipolarOutputSurroundRS,2))));

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));

% TODO: Calculate contrast gain adjustment

%% Apply rectification function and to the center and surround separately

% obj.rectify(input,'rType',{hw,fw,none})
obj.responseCenter   = obj.rectificationCenter(bipolarOutputLinearCenter);
obj.responseSurround = obj.rectificationSurround(bipolarOutputLinearSurround);

% Should we be rectifying the sum of the center/surround or should the two
% terms separately?
end