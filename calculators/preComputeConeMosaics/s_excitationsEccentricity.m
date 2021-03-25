%% Caclulate sum of excitations per unit area over different eccentricities
%
% Cone apertures increase, but the total number of cones decrease.  What is
% the next tradeoff as we measure at increasing eccentricities?
%
% BW/DHB
%
% See also
%  s_preComputeMosaics


%% Build a uniform scene
thisFOV = 9;

scene = sceneCreate('uniform');
scene = sceneSet(scene,'fov',thisFOV + 0.5);
oi = oiCreate; oi = oiCompute(oi,scene);

%% Count excitations as different eccentricities

fname = sprintf('cm-%ddeg-22-Mar-2021.mat',thisFOV);
fname = fullfile(isetRootPath,'local',fname);
foo = load(fname);
cm = foo.cm;
cm.integrationTime = 50e-3;
% cm.visualize;


%{
vParams = cm.visualize('params');
vParams.activation = noiseFree;
vParams.activationColorMap = hot(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);
%}

%%

c = 0:.5:(thisFOV-0.1);
centers = zeros(numel(c),2);
for ii=1:numel(c)
    centers(ii,:) = [c(ii) c(ii)]/2;
end

noiseFree = cm.compute(oi);

%%

fprintf('\n\n**********\n\n');
eccResponse = zeros(numel(c),1);
ecc = c(:);
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
    eccResponse(ii) = sum(roiResponse(:));
    fprintf('Center %.2f Excitations %.2f\n',centers(ii,1),sum(roiResponse(:)));
    
end

%% 
ieNewGraphWin;
plot(ecc,eccResponse,'-o','LineWidth',2);
xlabel('Eccentricity (deg)');
ylabel('Sum of excitations');
grid on

%%