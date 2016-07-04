function [spatialRFcenter, spatialRFsurround, rfDiaMagnitude, cellCenterLocations, tonicDrive] = buildSpatialRFArray(spacing, row, col, receptiveFieldDiameter1STDmicrons)
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
if spacing < 1e-2, spacing = spacing*1e6; end

patchSizeX = spacing;           % um
patchSizeY = spacing;           %(row/col)*spacing;
umPerSensorPx = patchSizeX/col; % CHANGE TO COL

%% Determine the number of RGCs in the mosaic

numberRGCsX = floor (patchSizeX / receptiveFieldDiameter1STDmicrons);
numberRGCsY = floor (patchSizeY / receptiveFieldDiameter1STDmicrons);

% Convert rf diameter in microns to rf diameter in pixels or cones
receptiveFieldDiameter1STD = receptiveFieldDiameter1STDmicrons/umPerSensorPx;

%% Set parameters for spatial RF difference of gaussians

% Modeling the RF shape as an ellipse, I guess?

% when extent gets to 5, the sum of the RF is ~0
extent = 2.5;

% Ellipse parameter
d1 = 1; d2 = 0;

% This is the quadratic form?
% d1 = 1; d2 = 0.25*randn(1,1);
Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

% Center has radius = 0.75*(radius of surround)
% r = .75; k = r;
r = 0.75; k = 1.032*r;

% Calculate coordiantes of 1 STD of center RF
% xv = rand(1,2);
% xvn = receptiveFieldDiameter1STD*xv./norm(xv);
% x1 = xvn(1); y1 = xvn(2);

% Need to do this for each RF? How to handle spatial asymmetry?
% magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);

%% Specify locations of RFs
rfctr = 0;

% Centers of RFs
icarr = 0+(0*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:(numberRGCsX-1)*receptiveFieldDiameter1STD);
jcarr = 0+(0*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:(numberRGCsY-1)*receptiveFieldDiameter1STD);

% A vector of all points out to the extent of the spatial RF
pts = (-extent*receptiveFieldDiameter1STD+1:1:extent*receptiveFieldDiameter1STD);

%% Create spatial RFs for each cell
tic
centerNoise = 1.25; % divide by 2 for mean offset

% Comment, please
centerCorrectY = 0+( 0+(jcarr(end) + 0 - (jcarr(1) + 0))/2 );  % + extent*receptiveFieldDiameter1STD;
centerCorrectX = 0+(0+ (icarr(end) + 0 - (icarr(1) + 0))/2 );  % + extent*receptiveFieldDiameter1STD;

% figure;   % What is this figure for?

iN = length(icarr); jN = length(jcarr);
cellCenterLocations = cell(iN,jN);
spatialRFArray    = cell(iN,jN);
spatialRFcenter   = cell(iN,jN);
spatialRFsurround = cell(iN,jN);
spatialRFonedim   = cell(iN,jN);  % One dimension
spatialRFFill     = cell(iN,jN);
rfDiaMagnitude    = cell(iN,jN);
spatialContours   = cell(iN,jN,2);
tonicDrive        = cell(iN,jN);

for icind = 1:length(icarr)
    
    for jcind = 1:length(jcarr)
        
        % Specify centers, offset even rows for hexagonal packing
        ic = icarr(icind) + centerNoise - (mod(jcind,2)-.5)*receptiveFieldDiameter1STD;
        jc = jcarr(jcind) - centerNoise;
        rfctr = rfctr+1;
   
        % Add some noise to deviate from circularity
        d1 = 1; d2 = 0;      % 0.0675*randn(1,1);
        Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);
        % receptiveFieldDiameter1STD == 1/sqrt(norm(Q)); % just to check

        % Calculate values for input to DoG function in an efficient way
        [i2, j2] = meshgrid(ic+pts,jc+pts);
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
        x1 = xvn(1); y1 = xvn(2);
        
        % magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        
        % Do some calculations to make plots where RFs are filled in
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{icind,jcind}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{icind,jcind,1} = magnitude1STD;
        
        % Get contours at 1STD
        cc = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
        % ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,1} = cc;
        
        clear cc
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr 
        cc = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
        % ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,2} = cc;
        
        rfDiaMagnitude{icind,jcind,2} = magnitude1STD;
        
        tonicDrive{icind,jcind} = 2.2702; % from ON Parasol 2013_08_19_6
    end
end
toc

end
