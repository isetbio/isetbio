function varargout = v_sceneoiSpatialResample(varargin)
% Spatial resample of a scene and an oi
%
% Validate the spatial resampling of the scene and oi photon data.
%
% BW, Copyright ISETBIO Team 2016
varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

function ValidationFunction(runTimeParams)
% Need to create validation parameters
%

ieInit;

%% Create a scene
scene = sceneCreate; 
scene = sceneSet(scene,'fov',1);

%% Resample
dx = 50;
scene = sceneSpatialResample(scene,dx,'um');
sr = sceneGet(scene,'spatial resolution','um');
assert(abs(sr(2) - 50) < 0.01)

if (runTimeParams.generatePlots)
    r = round(sceneGet(scene,'rows')/2);
    scenePlot(scene,'luminance hline',[1 r]); 
end
UnitTest.validationData('scene50', scene);


%% Create an oi with a larger spatial sample distance
scene = sceneCreate; 
scene = sceneSet(scene,'fov',10);
oi = oiCreate;
oi = oiCompute(oi,scene);
% ieAddObject(oi); oiWindow;

%% Now sample at 10 um
dx = 10;
oi = oiSpatialResample(oi,dx,'um');
sr = oiGet(oi,'spatial resolution','um');
assert(abs(sr(2) - dx) < 0.01)
if (runTimeParams.generatePlots)
    r = round(oiGet(oi,'rows')/2);
    oiPlot(oi,'illuminance hline',[1 r]); 
end
UnitTest.validationData('oi10', oi);

end
