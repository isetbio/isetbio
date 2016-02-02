function [spatialRFcenter, spatialRFsurround, rfDiaMagnitude, cellCenterLocations] = buildSpatialRFArray(spacing, row, col, receptiveFieldDiameter1STD)
% buildSpatialRF: a util function of the @rgc parent class that builds the
% spatial RF center and surround array.
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


patchSizeX = spacing; %sensorGet(sensor, 'width', 'um');
% sceneRows = sceneGet(scene,'rows');
sensorRows = row;% sensorGet(sensor,'rows');
umPerSensorPx = patchSizeX/sensorRows;

% % % NEED TO MAKE THIS DIFFERENT FOR OSLINEAR INITIALIZIATION!!!!

% coneSize = sensorGet(sensor, 'pixel size', 'um' );
patchSizeX = spacing;%sensorGet(sensor, 'width', 'um');
patchSizeY = spacing;%sensorGet(sensor, 'height', 'um');
% patchSizeX = sceneGet(scene,'sample size','um');
% patchSizeY = sceneGet(scene,'sample size','um');
% fov = sensorGet(sensor,'fov');
% numCones = sensorGet(sensor, 'size');


% Determine the number of RGCs in the mosaic
% numberRGCsX = floor (patchSizeX / receptiveFieldDiameter1STD);
% numberRGCsY = floor (patchSizeY / receptiveFieldDiameter1STD);
numberRGCsX = floor ((patchSizeX/umPerSensorPx) / receptiveFieldDiameter1STD);
numberRGCsY = floor ((patchSizeY/umPerSensorPx) / receptiveFieldDiameter1STD);

% Determine the RGC spacing in terms of number of cones.
% numberConesPerRF = floor (receptiveFieldDiameter1STD / coneSize(1));

% whene extent gets to 5, the sum of the RF is ~0
extent = 1.5;
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
for icind = 1:length(icarr)
    
    for jcind = 1:length(jcarr)
        ic = icarr(icind) - (mod(jcind,2)-1)*receptiveFieldDiameter1STD;% + centerNoise*(2*rand(1,1)-1);
        jc = jcarr(jcind);% + centerNoise*(2*rand(1,1)-1);
        rfctr = rfctr+1;
   
        d1 = 1; d2 = 0.0675*randn(1,1);
        Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);
        % receptiveFieldDiameter1STD == 1/sqrt(norm(Q)); % just to check

        [i2 j2] = meshgrid(ic+pts,jc+pts);
        i = i2(:); j = j2(:);
        
        IJ = bsxfun(@minus,[i j],[ic jc]);
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';
        %  icrm = repmat([ic jc],length(i),1);
        
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        so_center = reshape(exp(-0.5*p1),size(i2)); so_surround = reshape(k*exp(-0.5*p2),size(i2));
        so = so_center - so_surround;
        
        sx_cent = exp(-0.5*Q(1,1)*(0+pts).^2); sy_cent = exp(-0.5*Q(2,2)*(0+pts).^2);
        sx_surr = sqrt(k)*exp(-0.5*Q(1,1)*r*(0+pts).^2); sy_surr = sqrt(k)*exp(-0.5*Q(2,2)*r*(0+pts).^2);       
        
        cellCenterLocations{icind,jcind} = [ic jc];
%         load('rgc Parameters/OFFPar_1471_s1.mat');
%         load('rgc Parameters/OFFPar_1471_s2.mat');
%         so_center = abs(imresize(OFFPar_1471_s1,[length(pts),length(pts)]));
%         so_surround = abs(imresize(OFFPar_1471_s2,[length(pts),length(pts)]));
        spatialRFArray{icind,jcind} = so;
        spatialRFcenter{icind,jcind} = so_center;
        spatialRFsurround{icind,jcind} = so_surround;        

        spatialRFonedim{icind,jcind} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
                
        xv = rand(1,2);
        xvn = 1*receptiveFieldDiameter1STD*xv./norm(xv);
        % xvn = (extent-.1)*receptiveFieldDiameter1STD*xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
        
        % magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]);% - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{icind,jcind}  = find(abs(so_center)>magnitude1STD);
        rfDiaMagnitude{icind,jcind,1} = magnitude1STD;
        
        hold on;
        [cc,h] = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
%         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,1} = cc;
        
        clear cc h
        magnitude1STD = k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        % NOT SURE IF THIS IS RIGHT, bc contours are the same if so_surr 
        [cc,h] = contour(i2,j2,abs(so_center),[magnitude1STD magnitude1STD]);% close;
%         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind,2} = cc;
        
        rfDiaMagnitude{icind,jcind,2} = magnitude1STD;
    end
end
toc
close;
