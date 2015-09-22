function [spatialRFArray, spatialRFonedim, spatialContours, spatialRFFill, cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameter1STD)
% buildSpatialRF: a util function of the @rgc parent class
% 
% 
% 
% 
% 
% 
% 


coneSize = sensorGet(sensor, 'pixel size', 'um' );
patchSizeX = sensorGet(sensor, 'width', 'um');
patchSizeY = sensorGet(sensor, 'height', 'um');
fov = sensorGet(sensor,'fov');
numCones = sensorGet(sensor, 'size');


% Determine the number of RGCs in the mosaic
numberRGCsX = floor (patchSizeX / receptiveFieldDiameter1STD);
numberRGCsY = floor (patchSizeY / receptiveFieldDiameter1STD);

% Determine the RGC spacing in terms of number of cones.
numberConesPerRF = floor (receptiveFieldDiameter1STD / coneSize(1));


extent =3;
d1 = 1; d2 = 0;
% d1 = 1; d2 = 0.25*randn(1,1);
Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);

% kOn =  pi*325*0.03^2;
% kOff = pi*4.4*0.18^2;

r = .75; k = r;

xv = rand(1,2);
xvn = receptiveFieldDiameter1STD*xv./norm(xv);
x1 = xvn(1); y1 = xvn(2);
% Need to do this for each RF? How to handle spatial asymmetry?
magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);

% figure; 

tic
rfctr = 0;
icarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsX*receptiveFieldDiameter1STD;
jcarr = extent*receptiveFieldDiameter1STD:2*receptiveFieldDiameter1STD:numberRGCsY*receptiveFieldDiameter1STD;

% icarr = 1:2*receptiveFieldDiameter1STD:numberRGCsX*receptiveFieldDiameter1STD;
% jcarr = 1:2*receptiveFieldDiameter1STD:numberRGCsY*receptiveFieldDiameter1STD;


% s = zeros(ceil(max(icarr)+receptiveFieldDiameterParasol1STD),ceil(max(jcarr)+receptiveFieldDiameterParasol1STD));

% s = zeros(500);

pts = (-extent*receptiveFieldDiameter1STD+1:1:extent*receptiveFieldDiameter1STD);

tic

for icind = 1:length(icarr)
    ic = icarr(icind);
    for jcind = 1:length(jcarr)
        jc = jcarr(jcind);
        rfctr = rfctr+1;
   
        d1 = 1; d2 = 0.25*randn(1,1);
        Q = (1/receptiveFieldDiameter1STD^2)*[d1 d2; d2 d1]./norm([d1 d2; d2 d1]);


        [i2 j2] = meshgrid(ic+pts,jc+pts);
        i = i2(:); j = j2(:);
        
        IJ = bsxfun(@minus,[i j],[ic jc]);
        QIJ = Q*IJ'; rQIJ = r*Q*IJ';
        %  icrm = repmat([ic jc],length(i),1);
        
        p1 = prod([IJ(:,1) QIJ(1,:)'],2)+prod([IJ(:,2) QIJ(2,:)'],2);
        p2 = prod([IJ(:,1) rQIJ(1,:)'],2)+prod([IJ(:,2) rQIJ(2,:)'],2);
        so_center = exp(-0.5*p1); so_surround = k*exp(-0.5*p2);
        so = reshape(so_center - so_surround,size(i2));

        % s(floor(ic+(-extent*receptiveFieldDiameter1STD+1:extent*receptiveFieldDiameter1STD)),...
        %     floor(jc+(-extent*receptiveFieldDiameter1STD+1:extent*receptiveFieldDiameter1STD))) = so;
        
        sx_cent = exp(-0.5*Q(1,1)*(pts)); sy_cent = exp(-0.5*Q(2,2)*(pts));
        sx_surr = exp(-0.5*Q(1,1)*r*(pts)); sy_surr = exp(-0.5*Q(2,2)*r*(pts));       
        
        cellCenterLocations{icind,jcind} = [ic jc];
        spatialRFArray{icind,jcind} = so;
        spatialRFonedim{icind,jcind} = [(sx_cent - sx_surr); (sy_cent - sy_surr)];
                
        xv = rand(1,2);
        xvn = receptiveFieldDiameter1STD*xv./norm(xv);
        x1 = xvn(1); y1 = xvn(2);
        
%         magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
%         spatialRFFill{icind,jcind}  = find(so>magnitude1STD); 
        
        magnitude1STD = exp(-0.5*[x1 y1]*Q*[x1; y1]) - k*exp(-0.5*[x1 y1]*r*Q*[x1; y1]);
        spatialRFFill{icind,jcind}  = find(sx_cent>magnitude1STD);
        
        hold on;
        [cc,h] = contour(i2,j2,so,[magnitude1STD magnitude1STD]);% close;
%         ccCell{rfctr} = cc(:,2:end);
        cc(:,1) = [NaN; NaN];
        spatialContours{icind,jcind} = cc;
        
    end
end
toc
% figure; imagesc(sall(:,:,1)); shading flat
close;
% axis square
%%
% figure;
% for rfind = 1:rfctr
%     hold on; plot(ccCell{rfind}(1,2:end),ccCell{rfind}(2,2:end),'b');
% end
% title(sprintf('%s',namesCellTypes));