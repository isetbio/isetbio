function [sRFcenter, sRFsurround, rfDiaMagnitude, cellCenterLocations, tonicDrive] = ...
    buildSpatialRFArray(patchSize, inRow, inCol, rfDiameter)
% Builds the spatial RF center and surround for the cells in 
% 
% The spatial RFs are generated according to the size of the pixel, cone or
% bipolar mosaic, their spacing (in microns) and the diameter of the RGC RF
% as determined by the TEE of the retinal patch.
% 
%  [sRFcenter, sRFsurround, rfDiaMagnitude, cellCenterLocations] = 
%      buildSpatialRFArray(spacing, row, col, rfDiameter)
% 
% Inputs: 
%       spacing - Center to center of the RF in microns
%       row     - Number of input samples
%       col     - Number of input samples
%       rfDiameter - receptive field of 1 std in microns
%   
% Outputs: 
%       spatialRFcenter cell array, 
%       spatialRFsurround cell array,
%       rfDiaMagnitude at 1 std, 
%       cellCenterLocations cell array.
%  
%    WHAT IS TONIC DRIVE, and why is it here?
%         Tonic drive is the bias or DC term for the linear output of the
%         GLM. If the tonic drive term is greater than 0, then there is a
%         baseline firing rate even when the stimulus input is zero. It was
%         originally in the model and should not have been removed.
% 
% Example:
%   Build spatial RFs of the RGCs in this mosaic
%
% [obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
%     buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, obj.rfDiameter);
% 
% 9/2015 JRG (c) isetbio
% 7/2016 JRG updated

%% Manage parameters

% I think JRG built this assuming the patch size referred to a region of
% the cone mosaic.  
%
% Size of bipolar patch arrives in meters, which apparently is the width of
% the patch.  We convert to microns and build up a height/width
% representation (row,col)
patchSize = patchSize * 1e6;                      % Column (width) in um
patchSize = [patchSize*(inRow/inCol), patchSize]; % Row/Col in um

% Determine the number of RGC samples in the mosaic
% patchSize: um; rfDiameter: um/RGC cell; nRGC: RGC cells 
nRGC = floor(patchSize ./ rfDiameter); % number of rgc in h, v direction

% The rfDiameter comes here in units of um, and this 
% converts the rf diameter to units of number of bipolars
% rfDiameter in: um/RGC, 
% patchSize: um, inCol: number bipolar cells
% (patchSize(2) / inCol) : um/bipolar cell
% rfDiameter out: number bipolar cells per RGC
rfDiameter = rfDiameter / (patchSize(2) / inCol);

extent = 2.5;    % ratio between sampling size and spatial RF standard dev 
r = 0.75;        % radius ratio between center and surround for DoG
k = 1.032 * r;   % scaling of magnitude of surround

% centers of receiptive fields
centerNoise = 1.25;                                 % in number of bipolars; divide by 2 for mean offset
centerX = (0:2:nRGC(1)-1)*rfDiameter + centerNoise; % RGC center row coords in um
centerY = (0:2:nRGC(2)-1)*rfDiameter - centerNoise; % RGC center col coords in um
rows = length(centerX); cols = length(centerY);     % number of RGCs

% number bipolar cells out to the extent of the spatial RF
pts = -extent*rfDiameter+1 : extent*rfDiameter;

%% Create spatial RFs
% offset for zero centering in number of bipolars
centerCorrectY = (centerY(end) - (centerY(1)))/2; % um
centerCorrectX = (centerX(end) - (centerX(1)))/2; % um

% pre-allocate space
cellCenterLocations = cell(rows, cols);
spatialRFArray      = cell(rows, cols);
sRFcenter           = cell(rows, cols);
sRFsurround         = cell(rows, cols);
spatialRFonedim     = cell(rows, cols);
spatialRFFill       = cell(rows, cols);
rfDiaMagnitude      = cell(rows, cols);
% spatialContours     = cell(rows, cols, 2);
tonicDrive          = cell(rows, cols);

% s_center   =   exp(-0.5*(x-c)*Q*(x-c)') 
% s_surround = k*exp(-0.5*(x-c)*Q*(x-c)') 

% create spatial RFs for each cell
for ii = 1 : length(centerX)
    for jj = 1 : length(centerY)
        
        %% Compute 2D spatial RF
        % Specify centers in um, offset even rows for hexagonal packing
        ic = centerX(ii) - (mod(jj, 2) - 0.5) * rfDiameter + 3*centerNoise*(2*rand(1,1)-1);
        jc = centerY(jj) + 3*centerNoise*(2*rand(1,1)-1);
   
        % Add some noise to deviate from circularity 
        % (unitless: Q = (1/d^2)*[1 0; 0 1] yields circular SD with r = d
        d1 = 1; d2 =  10*0.0675*(rand(1,1)-0.5);      % 0.0675*randn(1,1);
        Q = (1/rfDiameter^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

        % Calculate values for input to DoG function in an efficient way
        [i2, j2] = meshgrid(ic+pts, jc+pts); % um
        i = i2(:); j = j2(:);                % um
        
        IJ = bsxfun(@minus,[i j],[ic jc]); % um
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';       % unitless
        %  icrm = repmat([ic jc],length(i),1);
        
        % (-0.5*(x-c)*Q*(x-c)'): unitless
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        
        % DoG calculation
        % conditional intensity, related by Poisson firing to spikes/sec
        so_center = reshape(exp(-0.5*p1), size(i2));
        so_surround = reshape(k*exp(-0.5*p2), size(i2));
        so = so_center - so_surround;
        
        %% Vectors (1D) instead of matrices (2D)        
        % conditional intensity, related by Poisson firing to spikes/sec
        sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2);
        sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
        sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2);
        sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);       
        
        % Store calculated parameters, units of conditional intensity
        cellCenterLocations{ii,jj} = [ic jc] - [centerCorrectX centerCorrectY]; % um
        spatialRFArray{ii,jj} = so;
        sRFcenter{ii,jj} = so_center;
        sRFsurround{ii,jj} = so_surround;        

        % units of conditional intensity
        spatialRFonedim{ii,jj} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
        
        %% Measure contour at 1 SD
        % Choose random orientation of asymmetrical RF to measure 1 SD
        % magnitude
        xv = rand(1,2);
        xvn = rfDiameter * xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
               
        % Do some calculations to make plots where RFs are filled in
        % Measure magnitude at 1 SD from center
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        % Find components of RF over that magnitude and store
        spatialRFFill{ii,jj}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{ii,jj,1} = magnitude1STD;        
        
        % clear cc
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);       
        
        rfDiaMagnitude{ii,jj,2} = magnitude1STD;
        
        %% Set tonic drive
        % Tonic drive is the bias or DC term for the linear output of the
        % GLM. If the tonic drive term is greater than 0, then there is a
        % baseline firing rate even when the stimulus input is zero.
        % Units of conditional intensity
        tonicDrive{ii,jj} = 2.2702; % from ON Parasol 2013_08_19_6
    end
end

end