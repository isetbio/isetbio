function [spatialRFcenter, spatialRFsurround, rfDiaMagnitude, cellCenterLocations] = buildSpatialRFArrayData(scene, sensor, receptiveFieldDiameter1STD)
% buildSpatialRFArrayData: a util function of the @rgc parent class that
% builds the spatial RF center and surround array with data from physiology
% collected in the Chichilnisky Lab.
% 
% Inputs: scene, sensor, and RF diameter.
%   
% Outputs: spatialRFcenter cell array, spatialRFsurround cell array,
% rfDiaMagnitude at 1 std, cellCenterLocations cell array.
% 
% Example:
% 
% % Build spatial RFs of all RGCs in this mosaic
% [obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
%     buildSpatialRFArray(scene, sensor, obj.rfDiameter);
% 
% 
% 
% (c) isetbio
% 9/2015 JRG
% 
% 

load('/Volumes/Lab/Users/james/isetbio/isettools/rgc/util/ej/data/2012-08-09-3/datarun.mat');
% Load only indices for on parasol
indices = [7 10 11 24 26 32 33 48 49 54 56 59 62 64 66 67 71 72 80  88 95 96 104 107 113 116 120 126 129 136 139 149 150 152 155 161 168 169 170 171 176 177 181 186 194 205 218 227 228 230 231 238 245 252];

patchSizeX = sensorGet(sensor, 'width', 'um');
sceneRows = sceneGet(scene,'rows');
umPerScenePx = patchSizeX/sceneRows;

% % % NEED TO MAKE THIS DIFFERENT FOR OSLINEAR INITIALIZIATION!!!!

coneSize = sensorGet(sensor, 'pixel size', 'um' );
patchSizeX = sensorGet(sensor, 'width', 'um');
patchSizeY = sensorGet(sensor, 'height', 'um');
% patchSizeX = sceneGet(scene,'sample size','um');
% patchSizeY = sceneGet(scene,'sample size','um');
% fov = sensorGet(sensor,'fov');
% numCones = sensorGet(sensor, 'size');


% Determine the number of RGCs in the mosaic
% numberRGCsX = floor (patchSizeX / receptiveFieldDiameter1STD);
% numberRGCsY = floor (patchSizeY / receptiveFieldDiameter1STD);
numberRGCsX = floor ((patchSizeX/umPerScenePx) / receptiveFieldDiameter1STD);
numberRGCsY = floor ((patchSizeY/umPerScenePx) / receptiveFieldDiameter1STD);

% Determine the RGC spacing in terms of number of cones.
% numberConesPerRF = floor (receptiveFieldDiameter1STD / coneSize(1));

% whene extent gets to 5, the sum of the RF is ~0
extent = 2;
d1 = 1; d2 = 0;
% d1 = 1; d2 = 0.25*randn(1,1);
Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

r = .75; k = r;

xv = rand(1,2);
xvn = receptiveFieldDiameter1STD*xv./norm(xv);
x1 = xvn(1); y1 = xvn(2);
% Need to do this for each RF? How to handle spatial asymmetry?
magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);

% figure; 

tic
rfctr = 0;

% icarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsX*receptiveFieldDiameter1STD;
% jcarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsY*receptiveFieldDiameter1STD;

icarr = 1*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsX*receptiveFieldDiameter1STD;
jcarr = 1*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsY*receptiveFieldDiameter1STD;

pts = (-extent*receptiveFieldDiameter1STD+1:1:extent*receptiveFieldDiameter1STD);

tic
centerNoise = 1.25; % divide by 2 for mean offset
for icind = 1:length(indices)
    
%     for jcind = 1:length(jcarr)
        ic = datarun.vision.sta_fits{indices(icind)}.mean(1); % icarr(icind);% + centerNoise*(2*rand(1,1)-1);
        jc = datarun.vision.sta_fits{indices(icind)}.mean(2); % jcarr(jcind);% + centerNoise*(2*rand(1,1)-1);
        
        rfctr = rfctr+1;
   
        d1 = datarun.vision.sta_fits{indices(icind)}.sd(2);
        d2 = datarun.vision.sta_fits{indices(icind)}.sd(1);
  %      Q = [d1 0; 0 d2];%./norm([d1 0; 0 d2]);
        Q = [d1 d2; d2 d1]./norm([d1 d2; d2 d1]);
        % Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);
        % receptiveFieldDiameter1STD == 1/sqrt(norm(Q)); % just to check

        
        %adjust orientation
        rfangle = -datarun.vision.sta_fits{indices(icind)}.angle;
        % rotate by angle and stretch
        RM = 1;%rmatrix2(rfangle / (2*pi) * 360);
        
        [i2 j2] = meshgrid(ic+pts,jc+pts);
        i = i2(:); j = j2(:);
        
        IJ = bsxfun(@minus,[i j],[ic jc]);
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';
        %  icrm = repmat([ic jc],length(i),1);
        
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        so_center = RM*reshape(exp(-0.5*(p1)),size(i2)); so_surround = RM*reshape(k*exp(-0.5*(p2)),size(i2));
        so = so_center - so_surround;
        
        sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2); sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
        sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2); sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);       
        
        cellCenterLocations{icind} = [ic jc];
        cx(icind) = ic; cy(icind) = jc;
%         load('rgc Parameters/OFFPar_1471_s1.mat');
%         load('rgc Parameters/OFFPar_1471_s2.mat');
%         so_center = abs(imresize(OFFPar_1471_s1,[length(pts),length(pts)]));
%         so_surround = abs(imresize(OFFPar_1471_s2,[length(pts),length(pts)]));
        spatialRFArray{icind} = so;
        spatialRFcenter{icind} = so_center;
        spatialRFsurround{icind} = so_surround;        

        spatialRFonedim{icind} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
                
        xv = rand(1,2);
        xvn = 1*1*xv./norm(xv);
        % xvn = (extent-.1)*receptiveFieldDiameter1STD*xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
        x1 = 1; y1 = 0;
        % magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{icind}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{icind,1} = magnitude1STD;
        
        hold on;
        [cc,h] = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
        % ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,1} = cc;
        
        clear cc h
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr 
        [cc,h] = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
         %ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,2} = cc;
        
        rfDiaMagnitude{icind,2} = magnitude1STD;
%     end
end
toc
close;
