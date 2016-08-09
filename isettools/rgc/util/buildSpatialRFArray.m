function [sRFcenter, sRFsurround, rfDiaMagnitude, cellCenterLocations, tonicDrive] = buildSpatialRFArray(spacing, inRow, inCol, rfDiameter)
% Builds the spatial RF center and surround arrays for each cell.
% 
% The spatial RFs are generated according to the size of the pixel, cone or
% bipolar mosaic, their spacing (in microns) and the diameter of the RGC RF
% as determined by the TEE of the retial patch.
% 
%   [sRFcenter, sRFsurround, rfDiaMagnitude, cellCenterLocations] = 
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
% Example:
% Build spatial RFs of all RGCs in this mosaic
% [obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
%     buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, obj.rfDiameter);
% 
% 9/2015 JRG (c) isetbio
% 7/2016 JRG updated

%% Manage parameters

% Spacing must be in microns
if spacing < 1e-2, spacing = spacing * 1e6; end
patchSize = [spacing*inRow/inCol spacing];  % width / height in um

% Determine the number of RGCs in the mosaic
nRGC = floor(patchSize ./ rfDiameter); % number of rgc in h, v direction

% Convert rf diameter in units of number of cones
% Notice this is based only columns, assuming the 
rfDiameter = rfDiameter / (patchSize(2) / inCol);

extent = 2.5;    % ratio between sampling size and spatial RF
r = 0.75;        % radius ratio between center and surround
k = 1.032 * r;   % 

% centers of receiptive fields
centerNoise = 1.25; % divide by 2 for mean offset
centerX = (0:2:nRGC(1)-1)*rfDiameter + centerNoise; 
centerY = (0:2:nRGC(2)-1)*rfDiameter - centerNoise;
rows = length(centerX); cols = length(centerY);

% points out to the extent of the spatial RF
pts = -extent*rfDiameter+1 : extent*rfDiameter;

%% Create spatial RFs
% offset for zero centering
centerCorrectY = (centerY(end) - (centerY(1)))/2;
centerCorrectX = (centerX(end) - (centerX(1)))/2;

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

% create spatial RFs for each cell
for ii = 1 : length(centerX)
    for jj = 1 : length(centerY)
        % Specify centers, offset even rows for hexagonal packing
        ic = centerX(ii) - (mod(jj, 2) - 0.5) * rfDiameter;
        jc = centerY(jj);
   
        % Add some noise to deviate from circularity
        d1 = 1; d2 = 0;      % 0.0675*randn(1,1);
        Q = (1/rfDiameter^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

        % Calculate values for input to DoG function in an efficient way
        [i2, j2] = meshgrid(ic+pts, jc+pts);
        i = i2(:); j = j2(:);
        
        IJ = bsxfun(@minus,[i j],[ic jc]);
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';
        %  icrm = repmat([ic jc],length(i),1);
        
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        
        % DoG calculation
        so_center = reshape(exp(-0.5*p1), size(i2));
        so_surround = reshape(k*exp(-0.5*p2), size(i2));
        so = so_center - so_surround;
        
        % Vectors instead of matrices
        sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2);
        sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
        sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2);
        sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);       
        
        % Store calculated parameters
        cellCenterLocations{ii,jj} = [ic jc] - [centerCorrectX centerCorrectY];
        spatialRFArray{ii,jj} = so;
        sRFcenter{ii,jj} = so_center;
        sRFsurround{ii,jj} = so_surround;        

        spatialRFonedim{ii,jj} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
                
        xv = rand(1,2);
        xvn = rfDiameter * xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
               
        % Do some calculations to make plots where RFs are filled in
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{ii,jj}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{ii,jj,1} = magnitude1STD;        
        
        % clear cc
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);       
        
        rfDiaMagnitude{ii,jj,2} = magnitude1STD;
        
        tonicDrive{ii,jj} = 2.2702; % from ON Parasol 2013_08_19_6
    end
end

end