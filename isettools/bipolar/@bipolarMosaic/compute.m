function [obj, nTrialsCenter, nTrialsSurround] = compute(obj, cmosaic, varargin)
% BIPOLAR.COMPUTE - Compute bipolar continuous current responses
%
%    bipolar.compute(coneMosaic,varargin);
%
% The bipolars act as a spatial-temporal function that converts the cone
% photocurrent into bipolar current that is delivered to the retinal
% ganglion cells.
%
% This could be transformed to be like the rgc computation where each
% cell has its own spatial RF.  Not there yet, just using convolutions.
%
%
% Required parameters:
%   obj:       a bipolar object
%   cmosaic:   coneMosaic (current must be computed)
%
% Bipolar cell parameters are set when the bp mosaic is created.  Important
% parameters include:
%
%  'Cell type': The bipolar cells are classified into several types
%     * on/off diffuse, which connect to parasol RGCs
%     * on/off midget, which connect to midget RGCs
%     * on small bistratified, which connect to S-cone bistratified
%  'Processing' - linear or rectified
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
%  *** This feature has not been tested or much used - on our list
%  REFERENCE: Meister option
%
%  The spatial filtering is followed by a temporal filter that is selected
%  in order to match the impulse response that expect to find at the RGC
%  input.  REFERENCE: Chichilnisky option
%
% Returns
%
%  The bipolar object handles multiple trials, and these are returned
%  separately for the center and surround of the biplar model when use have
%  output arguments, as below.
%
%    [~, bpNTrialsCenter, bpNTrialsSurround] = ...
%           bp.compute(cMosaic,'nTrialsInput',cMosaicCurrentNTrials);
%
%  The cMosaicCurrentNTrials has dimensions (nTrials,row,col,nFrames), as
%  returned by cmosaic.computeCurrent.  The full bipolar response simply
%  add bpNTrialsCenter + bpNTrialsSurround.
%
%  Only the last trial reponse is stored in the object itself.
%
% 5/2016 JRG (c) isetbio team

%% parse input parameters

p = inputParser;
p.addRequired('obj', @(x) (isa(x, 'bipolarMosaic')));
p.addRequired('cmosaic', @(x) (isa(x, 'coneMosaic')));
addParameter(p, 'coneTrials',  [], @isnumeric);

% parse - no options at this opint
p.parse(obj, cmosaic, varargin{:});

coneTrials = p.Results.coneTrials;

if isempty(cmosaic.current)
    error('No cone photocurrent.  Use cmosaic.computeCurrent.');
end

if ~isempty(coneTrials),  nTrials = size(coneTrials,1);
else,                     nTrials = 1;
end

%% Spatial filtering and subsampling

% Convolve spatial RFs across the photo current of the cones in the mosaic

for iTrial = 1:nTrials
    
    % This places the cone 3D matrix into a coneNumber x time matrix
    if ~isempty(coneTrials)
        osSig = RGB2XWFormat(squeeze(coneTrials(iTrial,:,:,:)));
    else
        osSig = RGB2XWFormat(cmosaic.current);
    end
    
    
    %% Enfoce anatomical rules on cone to bipolar connections
    
    % Anatomical rules:
    %
    %  off Diffuse, on Diffuse and on Midget - These receive no S cone input
    %  offMidget - keep S cones but scale the connection strength down by 75%
    %  onSBC     - S cone inputs to center, only L/M cone inputs to surround
    %
    % Citations:  See bipolar.m.  Wiki page <>
    
    % TODO:  We should set this up as a structure that we use to implement
    % the anatomical rules.  Let's send in a struct that defines the
    % anatomical rules (e.g., anatomyRules) with slots that implement the
    % kind of stuff listed above.
    %
    switch obj.cellType
        case{'offdiffuse','ondiffuse','onmidget'}
            
            osSigCenter   = osSig;
            osSigSurround = osSig;
            
            % Remove S cone input for these types of bipolars
            
            % Find the locations indices of the different cone types
            [~,~,S] = coneTypeLocations(cmosaic,'format','index');
            
            % Zero the photocurrent of the S cones. Do this for both the center
            % and the surround.
            
            minval = min(osSig(:));
            osSigCenter(S(:),:)   = minval*ones(size(osSigCenter(S,:)));
            osSigSurround(S(:),:) = minval*ones(size(osSigCenter(S,:)));
            
        case{'offmidget'}
            % Keep S cone input for off Midget but only weight by 0.25
            
            % Find the locations (row, col) of the different cone types
            [~,~,S] = coneTypeLocations(cmosaic,'format','index');
            
            minval = min(osSig(:));
            
            osSigCenter   = osSig;
            osSigCenter(S,:)   = 0.25*(osSigCenter(S,:)-minval)+minval;
            
            osSigSurround   = osSig;
            osSigSurround(S,:) = 0.25*(osSigSurround(S,:)-minval)+minval;
            
        case{'onsbc'}
            % Set L and M cones to zero in SBC center, set S cones to zero
            % in SBC surround.
            
            % Find the indices of the different cone types
            [L,M,S] = coneTypeLocations(cmosaic,'format','index');
            
            % This is one long vector of L,M cone indices
            LM = [L; M];
            
            % Find the effectively zero outer segment signal for this
            % mosaic
            minval = min(osSig(:));
            
            % When the center is an LM cone, make all of the time steps in
            % the center the smallest value
            osSigCenter       = osSig;
            osSigCenter(LM,:) = minval*ones(size(osSigCenter(LM,:)));
            
            % Put effectively zero S-cone signals into the surround
            osSigSurround      = osSig;
            osSigSurround(S,:) = minval*ones(size(osSigSurround(S,:)));
            
    end
    
    % Put the data back into RGB format, like RGB2XW()
    sz = size(cmosaic.current);
    osSigCenter   = XW2RGBFormat(osSigCenter,sz(1),sz(2));
    osSigSurround = XW2RGBFormat(osSigSurround,sz(1),sz(2));

    % cmosaic.window;
    % vcNewGraphWin; ieMovie(osSigCenter);
    
    %% Spatial filtering 
    
    % Full spatial convolution for every frame.  The kernel is only 2D
    % which is why we have a space-only convolution.
    bipolarCenter   = ieSpaceTimeFilter(osSigCenter, obj.sRFcenter);
    bipolarSurround = ieSpaceTimeFilter(osSigSurround, obj.sRFsurround);
    % vcNewGraphWin; ieMovie(bipolarCenter);
    % vcNewGraphWin; ieMovie(bipolarSurround);

    % The bipolar cells might not be abutting, in some model.  We have
    % never used that condition - they always fully tile.  So spacing is
    % always 1 and we haven't subsampled.  This is here because, well, we
    % might some day. (BW/JRG).
    spacing = size(obj.sRFcenter,1);
    if spacing ~= 1
        % Subsample in space to the resolution for this bipolar mosaic. The
        % spacing is equal to the number of pixels that make up the center
        % of the spatial receptive field.  This could be a settable
        % parameter for others to experiment with, too.  We need a
        % reference.
        bipolarCenter   = ieImageSubsample(bipolarCenter, spacing);
        bipolarSurround = ieImageSubsample(bipolarSurround, spacing);
    end
    
    %% Temporal filtering
    
    % Reshape the data for the temporal convolution
    [bipolarCenter, row, col] = RGB2XWFormat(bipolarCenter);
    bipolarSurround = RGB2XWFormat(bipolarSurround);
    
    % This is the impulse response filter
    bipolarFilt = bipolarFilter(obj, cmosaic);
    
    % If we wanted to rectify the signal, we could do it here
    %
    % obj.rectify(input,'rType',{hw,fw,none})
    % obj.responseCenter   = obj.rectificationCenter(bipolarOutputLinearCenter);
    % obj.responseSurround = obj.rectificationSurround(bipolarOutputLinearSurround);
    %
    
    %% Rectification and temporal convolution issues
    
    % Rectification - not tested or analyzed 
    
    % We have in the past shifted the bipolar response levels to a minimum
    % of zero.  That is arbitrary and produces higher contrast signals.  It
    % might be OK because we have no real units on the bipolar current.
    % Or, maybe we should leave them alone.  Anyway, here we are shifting
    % them to a min of zero.
    
    % bipolarCenter = obj.rectificationCenter(bipolarCenter - (min(bipolarCenter')'*ones(1,size(bipolarCenter,2))));
    bipolarCenter = bipolarCenter - (min(bipolarCenter')'*ones(1,size(bipolarCenter,2)));
    tmpCenter = conv2(bipolarFilt,bipolarCenter);
    % vcNewGraphWin; tmp = XW2RGBFormat(tmpCenter,row, col); ieMovie(tmp);
    
    % Rectification
    % Not fully tested or analyzed - 
    % bipolarSurround = obj.rectificationSurround(bipolarSurround-(min(bipolarSurround')'*ones(1,size(bipolarSurround,2))));
    bipolarSurround = bipolarSurround-(min(bipolarSurround')'*ones(1,size(bipolarSurround,2)));
    tmpSurround = conv2(bipolarFilt,bipolarSurround);
    % vcNewGraphWin; tmp = XW2RGBFormat(tmpSurround,row, col); ieMovie(tmp);

    if ~isempty(coneTrials)
        if iTrial == 1
            nTrialsCenter = zeros([nTrials,size(XW2RGBFormat(tmpCenter(:,1:cmosaic.tSamples),row,col))]);
            nTrialsSurround = zeros([nTrials,size(XW2RGBFormat(tmpSurround(:,1:cmosaic.tSamples),row,col))]);
        end
        
        nTrialsCenter(iTrial,:,:,:) = XW2RGBFormat(tmpCenter(:,1:cmosaic.tSamples),row,col);
        nTrialsSurround(iTrial,:,:,:) = XW2RGBFormat(tmpSurround(:,1:cmosaic.tSamples),row,col);
        %     obj.responseCenter(iTrial,:,:,:) = tmpTrialCenter;
        %     obj.responseSurround(iTrial,:,:,:) = tmpTrialSurround;
        
    else
        nTrialsCenter   = 0;
        nTrialsSurround = 0;
    end
    
    if iTrial == nTrials
        obj.responseCenter   = XW2RGBFormat(tmpCenter(:,1:size(cmosaic.current,3)),row,col);
        obj.responseSurround = XW2RGBFormat(tmpSurround(:,1:size(cmosaic.current,3)),row,col);
    end
    
    
end

end