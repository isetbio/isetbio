function absBlurred = layerSpatialConv(layer)
% Spatial convolution of the cone voltage with each of the RGC components
%
%   absBlurred = layerSpatialConv(layer)
%
% Convolve the temporal samples of the cone volts (cVolts) with each
% component receptive field of this layer.  The center and surround are
% kept separate because they subsequently are processed by different time
% courses.
%
% absBlurred is a 4D matrix corresponding to the (x,y,t,component)
%
% The units of absBlurred are volts.
%
% Example:
%
% See also:  rgcComputeSpikes, layerTemporalConv
%
% (c) Stanford VISTA Team 2010

%% Check input
if notDefined('layer'), error('Input a rgcLayer object'); end

%% Extract values

rgcP        = layer.get('parent');
nFrames     = rgcP.get('nFrames');      % Image dimensions

% Image values
cVolts      = rgcP.get('cone voltages');  % Cone cVolts
sensor      = rgcP.get('sensor');
r = size(cVolts,1); c = size(cVolts,2);

% RF components are the receptive field profiles for the L,M and S cones
% for the different components, say center and surround.
% RF(r,c,component,coneType)
RF    = layer.get('RF components');   % Spatial profile
nComponent = size(RF,3);
coneWgts =  layer.get('cone weights'); % ones(nComponent,3);        % Always positive?
% coneWgts = layer.get('RF coeffs');  % Weights for each cone type

% Allocate space for blurred output.  Not yet downsampled.
bigBlurred = zeros(r,c,nFrames,nComponent);

% The cones are stored as K,L,M,S (Black and then the LMS). So the jj loop
% from 1:3 has to pull out the cone images with a plus one (see below).
for ii=1:nFrames
    lms = plane2rgb(cVolts(:,:,ii),sensor,0);
    for cc = 1:nComponent
        cv = zeros(r,c);
        for jj=1:3
            cv = cv + coneWgts(cc,jj)*conv2(lms(:,:,(jj+1)), RF(:, :, cc), 'same');
        end
        % vcNewGraphWin; imagesc(cv); colormap(gray)
        bigBlurred(:,:,ii,cc) = cv;    
    end
end
% tmp = squeeze(sum(bigBlurred,3));
% vcNewGraphWin;  mesh(tmp(:,:,1)); 
% vcNewGraphWin;  mesh(tmp(:,:,2)); 
% vcNewGraphWin;  mesh(tmp(:,:,1) + tmp(:,:,2)); 

% Next, do the spatial subsampling using cell and cone spacing.
% Should this be rounded?  We could also interpolate rather than sample.
samplePositions = layer.get('cellConePos');

% cellSpacing = round(layer.get('cell spacing') / rgcP.get('cone spacing'));
% rS = 1:cellSpacing:r;
% cS = 1:cellSpacing:c;

absBlurred = bigBlurred(samplePositions{1},samplePositions{2},:,:);

% tmp = squeeze(sum(absBlurred,3));
% vcNewGraphWin;  mesh(tmp(:,:,1)); 
% vcNewGraphWin;  mesh(tmp(:,:,2)); 
% vcNewGraphWin;  mesh(tmp(:,:,1) + tmp(:,:,2)); 

% Receptive fields define weighted sums of the input voltages from the
% cones.  The first rf component in a center-surround  sums to 1 and thus
% preserves the mean.  The 2nd component (usually the surround) has a
% negative value and subtracts away some of the mean (and modifies the
% contrast).  So sum(RF(:)) in this case is < 1.

return

