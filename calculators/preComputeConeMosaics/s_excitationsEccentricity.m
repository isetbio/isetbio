%% Count excitations as different eccentricities

thisFOV = 1;
fname = sprintf('cm-%ddeg-22-Mar-2021.mat',thisFOV);
fname = fullfile(isetRootPath,'local',fname);
foo = load(fname);
cm = foo.cm;
cm.visualize;

vParams = cm.visualize('params');
cm.visualize(vParams);
cm.visualize;

scene = sceneCreate('uniform');
scene = sceneSet(scene,'fov',6);
oi = oiCreate; oi = oiCompute(oi,scene);
cm.integrationTime = 50e-3;
noiseFree = cm.compute(oi);

vParams = cm.visualize('params');
vParams.activation = noiseFree;
vParams.activationColorMap = hot(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);

%%

c = 0:.5:(thisFOV-1);
centers = [];
for ii=1:numel(c)
    centers(end+1,:) = [c(ii) c(ii)]/2;
end

%%

fprintf('\n\n**********\n\n');

for ii=1:size(centers,1)
    
    roi = struct('center', centers(ii,:), 'radius', 0.1);
    
    % Compute border of ROI
    roiBorderX = roi.center(1) + roi.radius*cosd(0:10:360);
    roiBorderY = roi.center(2) + roi.radius*sind(0:10:360);
    % ieNewGraphWin; plot(roiBorderX,roiBorderY)
    
    % Find indices of cones within the ROI border
    [in,on] = inpolygon(cm.coneRFpositionsDegs(:,1),cm.coneRFpositionsDegs(:,2),roiBorderX,roiBorderY);
    indicesOfConesWithinROI = in | on;
    
    % Form a response vector in which only cones within the ROI have a non-zero response
    roiResponse = 0*noiseFree;
    roiResponse(indicesOfConesWithinROI) = noiseFree(indicesOfConesWithinROI);
    fprintf('Center %.2f Excitations %.2f\n',centers(ii,1),sum(roiResponse(:)));
end

