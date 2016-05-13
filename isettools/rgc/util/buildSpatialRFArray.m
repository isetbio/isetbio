function [spatialRFcenter, spatialRFsurround, rfDiaMagnitude, cellCenterLocations] = buildSpatialRFArray(spacing, row, col, receptiveFieldDiameter1STDmicrons)
%% buildSpatialRF builds the spatial RF center and surround arrays for each cell
% The spatial RFs are generated according to the number of pixel or cone
% inputs, their spacing (in microns) and the diameter of the RF as
% determined by the TEE of the retial patch.
% 
%   [spatialRFcenter, spatialRFsurround, rfDiaMagnitude, cellCenterLocations] = 
%             buildSpatialRFArray(spacing, row, col, receptiveFieldDiameter1STDmicrons)
% 
%          %  [only called internally from @rgcMosaic/initalize.m]
% 
% Inputs: 
%       spacing, 
%       row, 
%       col, 
%       receptiveFieldDiameter1STDmicrons.
%   
% Outputs: 
%       spatialRFcenter cell array, 
%       spatialRFsurround cell array,
%       rfDiaMagnitude at 1 std, 
%       cellCenterLocations cell array.
% 
% Example:
% 
% See @rgcMosaic/initialize.m
% Build spatial RFs of all RGCs in this mosaic
% [obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
%     buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, obj.rfDiameter);
% 
% 9/2015 JRG (c) isetbio

%% Find number of pixels/cones per RGC spatial RF
% Calculate microns/pixels or microns/cone

% Check spacing is in meters
if spacing < 1e-2
    spacing = spacing*1e6;
end

% patchSizeX = spacing; % um 
% patchSizeY = spacing; % um 
% sensorRows = row;     % 
% umPerSensorPx = patchSizeX/sensorRows;

patchSizeX = spacing;
patchSizeY = spacing; %(row/col)*spacing;
% sensorRows = row;  
umPerSensorPx = patchSizeX/col; % CHANGE TO COL
%% Determine the number of RGCs in the mosaic
numberRGCsX = floor (patchSizeX / receptiveFieldDiameter1STDmicrons);
numberRGCsY = floor (patchSizeY / receptiveFieldDiameter1STDmicrons);
% numberRGCsX = floor ((patchSizeX/umPerSensorPx) / receptiveFieldDiameter1STD);
% numberRGCsY = floor ((patchSizeY/umPerSensorPx) / receptiveFieldDiameter1STD);

% Convert rf diameter in microns to rf diameter in pixels or cones
receptiveFieldDiameter1STD = receptiveFieldDiameter1STDmicrons/umPerSensorPx;

% Determine the RGC spacing in terms of number of cones.
% numberConesPerRF = floor (receptiveFieldDiameter1STD / coneSize(1));

%% Set parameters for spatial RF difference of gaussians
% whene extent gets to 5, the sum of the RF is ~0
extent = 2.5;

% Ellipse parameter
d1 = 1; d2 = 0;
% d1 = 1; d2 = 0.25*randn(1,1);
Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

% Center has radius = 0.75*(radius of surround)
r = .75; k = r;

% Calculate coordiantes of 1 STD of center RF
xv = rand(1,2);
xvn = receptiveFieldDiameter1STD*xv./norm(xv);
x1 = xvn(1); y1 = xvn(2);
% Need to do this for each RF? How to handle spatial asymmetry?
magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);

%% Specify locations of RFs
tic
rfctr = 0;

% icarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsX*receptiveFieldDiameter1STD;
% jcarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsY*receptiveFieldDiameter1STD;

% Centers of RFs
icarr = 0+(0*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:(numberRGCsX-1)*receptiveFieldDiameter1STD);
jcarr = 0+(0*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:(numberRGCsY-1)*receptiveFieldDiameter1STD);

% A vector of all points out to the extent of the spatial RF
pts = (-extent*receptiveFieldDiameter1STD+1:1:extent*receptiveFieldDiameter1STD);
% pts = (0:1:(extent*receptiveFieldDiameter1STD)); pts = .5+[-pts(end:-1:2) pts];

%% Create spatial RFs for each cell
tic
centerNoise = 0;%1.25; % divide by 2 for mean offset

% centerCorrectY = 0+( 0+(jcarr(end) + pts(end) - (jcarr(1) + pts(1)))/2 )% - receptiveFieldDiameter1STD/4;
% centerCorrectX = 0+(0+ (icarr(end) + pts(end) - (icarr(1) + pts(1)))/2 )

centerCorrectY = 0;%+( 0+(jcarr(end) + 0 - (jcarr(1) + 0))/2 );% + extent*receptiveFieldDiameter1STD;
centerCorrectX = 0;%+(0+ (icarr(end) + 0 - (icarr(1) + 0))/2 );% + extent*receptiveFieldDiameter1STD;


figure;

for icind = 1:length(icarr)
    
    for jcind = 1:length(jcarr)
        % Specify centers, offset even rows for hexagonal packing
        ic = icarr(icind) + centerNoise*(2*rand(1,1)-1) - (mod(jcind,2)-.5)*receptiveFieldDiameter1STD;
        jc = jcarr(jcind) + centerNoise*(2*rand(1,1)-1);
        rfctr = rfctr+1;
   
        % Add some noise to deviate from circularity
        d1 = 1; d2 = 0;%0.0675*randn(1,1);
        Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);
        % receptiveFieldDiameter1STD == 1/sqrt(norm(Q)); % just to check

        % Calculate values for input to DoG function in an efficient way
        [i2 j2] = meshgrid(ic+pts,jc+pts);
        i = i2(:); j = j2(:);
        
        IJ = bsxfun(@minus,[i j],[ic jc]);
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';
        %  icrm = repmat([ic jc],length(i),1);
        
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        
        % DoG calculation
        so_center = reshape(exp(-0.5*p1),size(i2)); so_surround = reshape(k*exp(-0.5*p2),size(i2));
        so = so_center - so_surround;
        
        % Vectors instead of matrices
        sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2); sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
        sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2); sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);       
        
        % Store calculated parameters
        cellCenterLocations{icind,jcind} = [ic jc] - [centerCorrectX centerCorrectY];
        
        spatialRFArray{icind,jcind} = so;
        spatialRFcenter{icind,jcind} = so_center;
        spatialRFsurround{icind,jcind} = so_surround;        

        spatialRFonedim{icind,jcind} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
                
        xv = rand(1,2);
        xvn = 1*receptiveFieldDiameter1STD*xv./norm(xv);
        % xvn = (extent-.1)*receptiveFieldDiameter1STD*xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
        
        % magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        
        % Do some calculations to make plots where RFs are filled in
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{icind,jcind}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{icind,jcind,1} = magnitude1STD;
        
        hold on;
        % Get contours at 1STD
        [cc,h] = contour(i2-0,j2-0,abs(so_center),[magnitude1STD magnitude1STD]);% close;
        % ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,1} = cc;
        
        clear cc h
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr 
        [cc,h] = contour(i2-0,j2-0,abs(so_center),[magnitude1STD magnitude1STD]);% close;
        % ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,2} = cc;
        
        rfDiaMagnitude{icind,jcind,2} = magnitude1STD;
    end
end
toc
close;
