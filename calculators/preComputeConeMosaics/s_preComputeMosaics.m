%% Pre-calculate an array of cone mosaics
%
%  Standard method on March 21 

% A number of field of views
fov = 1:10;

for ii=1:5
    cm = cMosaic('sizeDegs', [fov(ii),fov(ii)]);
    fname = sprintf('cm-%ddeg-%s.mat',fov(ii),date);
    fname = fullfile(isetRootPath,'local',fname)
    save(fname,'cm');
end

%%

thisFOV = 3;
fname = sprintf('cm-%ddeg-%s.mat',thisFOV,date);
fname = fullfile(isetRootPath,'local',fname);
foo = load(fname);
cm = foo.cm;

scene = sceneCreate('ringsrays');
scene = sceneSet(scene,'fov',5);
oi = oiCreate; oi = oiCompute(oi,scene);
cm.integrationTime = 50e-3;
noiseFree = cm.compute(oi);

vParams = cm.visualize('params');
vParams.activation = noiseFree;
vParams.activationColorMap = hot(512);
vParams.verticalActivationColorBar = true;
vParams.activationRange = [0 max(noiseFree(:))];
cm.visualize(vParams);
