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

%% TODO NOTES
%
% Load bipolar temporal filter

% Bipolar filters were deconvolved from the measured temporal impulse response of each
% cell in the mosaic and the linear cone temporal response. The
% mean of the bipolar temporal filters for the whole mosaic is used
% as the ideal bipolar filter.
% 
% There are different bipolar filters for on p/m and off p/m cells
% if strcmpi(obj.cellType,'offDiffuse')
%     % Off parasol (off diffuse) only
%     data = load([isetRootPath '/data/bipolar/bipolarFilt_200_OFFP_2013_08_19_6_all.mat']);
% else
%     % On parasol and all others
%     data = load([isetRootPath '/data/bipolar/bipolarFilt_200_ONP_2013_08_19_6_all.mat']);
% end
% The bipolarFiltMat is 100 different filters.  Fix that in the
% representation.  The script that builds this is about to be uploaded to
% the repository by JRG.
% bipolarFilt = -mean(data.bipolarFiltMat)';
% TODO: Remove bipolar.filterType parameter from object
%
%
% BW removed: Zero-mean the input signal.
%
% Why?
%
% Shouldn't the receptive field should govern how we map the input current
% to the output current (BW)?   If the RF has a zero mean, then spatial
% mean maps to 0. If the RF has a unit mean then spatial mean maps to mean
%
% if size(osSig,2) > 1
%     % Typical case.  Substract the mean over time of each cone signal from
%     % itself.
%     osSig = bsxfun(@minus, osSig, mean(osSig, 2));
% end


%% parse input parameters

p = inputParser;
p.addRequired('obj', @(x) (isa(x, 'bipolar')));
p.addRequired('cmosaic', @(x) (isa(x, 'coneMosaic')));  

% parse - no options at this opint
p.parse(obj, cmosaic, varargin{:});

if isempty(cmosaic.current), 
    error('No cone photocurrent.  Use cmosaic.computeCurrent.'); 
end
%% Spatial filtering and subsampling

% Convolve spatial RFs across the photo current of the cones in the mosaic

% This places the cone 3D matrix into a coneNumber x time matrix
osSig = RGB2XWFormat(cmosaic.current);


%% Enfoce anatomical rules on cone to bipolar connections

% Anatomical rules:
%
%  off Diffuse, on Diffuse and on Midget - These receive no S cone input
%  offMidget - keep S cones but scale the connection strength down by 75% 
%  onSBC     - S cone inputs to center, only L/M cone inputs to surround
%
% Citations:  See bipolar.m.  Wiki page <>

% TODO:  We should set this up as a structure that we use to implement the
% anatomical rules.  Let's send in a struct that defines the anatomical
% rules (e.g., aRules) with slots that implement the kind of stuff listed
% above.
%
switch obj.cellType
    case{'offdiffuse','ondiffuse','onmidget'}

        osSigCenter   = osSig;
        osSigSurround = osSig;

        % Remove S cone input for these types of bipolars

        % Find the locations (row, col) of the different cone types
        [~,~,S] = coneTypeLocations(cmosaic,'val','index');
        
        % Zero the photocurrent of the S cones. Do this for both the center
        % and the surround.
        z = zeros(length(S),size(osSig,2));
        osSigCenter(S(:),:)   = z;
        osSigSurround(S(:),:) = z;
        
    case{'offmidget'}
        % Keep S cone input for off Midget but only weight by 0.25
        
        % Find the locations (row, col) of the different cone types
        [~,~,S] = coneTypeLocations(cmosaic,'val','index');
        
        minval = min(osSig(:));
        
        osSigCenter   = osSig;
        osSigCenter(S,:)   = 0.25*(osSigCenter(S,:)-minval)+minval;
        
        osSigSurround   = osSig;
        osSigSurround(S,:) = 0.25*(osSigSurround(S,:)-minval)+minval;

    case{'onsbc'}  
        % Set L and M cones to zero in SBC center, set S cones to zero in
        % SBC surround.
        % Find the locations (row, col) of the different cone types
        [L,M,S] = coneTypeLocations(cmosaic,'val','index');
        LM = [L; M];
                
        minval = min(osSig(:));
        % Set center to only have S cones
        
        osSigCenter   = osSig;
        osSigCenter(LM,:)   = minval*ones(size(osSigCenter(LM,:)));
        
        osSigSurround   = osSig;
        osSigSurround(S,:)   = minval*ones(size(osSigSurround(S,:)));
                       
end

% Put the data back into RGB format, like RGB2XW()
osSigCenter   = reshape(osSigCenter,size(cmosaic.current));
osSigSurround = reshape(osSigSurround,size(cmosaic.current));

%% Spatial convolution

% Full spatial convolution for every frame
bipolarCenter   = ieSpaceTimeFilter(osSigCenter, obj.sRFcenter);
bipolarSurround = ieSpaceTimeFilter(osSigSurround, obj.sRFsurround);

% Subsample in space to the resolution for this bipolar mosaic.
% The spacing is equal to the number of pixels that make up the center of
% the spatial receptive field.  This could be a settable parameter for
% others to experiment with, too.  We need a reference.
spacing = size(obj.sRFcenter,1);
bipolarCenter   = ieImageSubsample(bipolarCenter, spacing);
bipolarSurround = ieImageSubsample(bipolarSurround, spacing);

%% Temporal filtering

% Reshape for temporal convolution
[bipolarCenter, row, col] = RGB2XWFormat(bipolarCenter);
bipolarSurround = RGB2XWFormat(bipolarSurround); 

%% New method

% The filter isn't right.  Time base is off.  Let's deal with it.
bipolarFilt = bipolarFilter(obj, cmosaic);

%% Compute the temporal response of the bipolar mosaic
%
% Deal with rectification, next. Need a plan
%
% obj.rectify(input,'rType',{hw,fw,none})
% obj.responseCenter   = obj.rectificationCenter(bipolarOutputLinearCenter);
% obj.responseSurround = obj.rectificationSurround(bipolarOutputLinearSurround);
%

tmp = conv2(bipolarFilt,bipolarCenter);
obj.responseCenter = XW2RGBFormat(tmp(:,1:cmosaic.tSamples),row,col);

tmp = conv2(bipolarFilt,bipolarSurround);
obj.responseSurround = XW2RGBFormat(tmp(:,1:cmosaic.tSamples),row,col);

end