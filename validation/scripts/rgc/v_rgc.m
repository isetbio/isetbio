function varargout = v_rgc(varargin)
%
% Validate the RGC object.


    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;



%% Compare spatial response to spatial receptive field and
%%% temporal impulse response to temporal response with impulse stimulus

%%% build scene
params.image_size = 63;
params.meanLuminance = 100;
params.nsteps = 1;
params.fov = 0.8;
[scene, display] = sceneHorwitzHassWhiteNoise(params);

% sceneRGB = sceneHorwitzHassWhiteNoiseRGB(params);
% sceneRGB = sceneHorwitzHassBarRGB(params);
sceneRGB = zeros(params.image_size, params.image_size, params.nsteps, 3);
sceneRGB(23,23,1,:) = [1 1 1];
% sceneRGB(38,38,1,:) = [1 1 1];

oi  = oiCreate('wvf human');
%%% build sensor for white noise
sensor = sensorHorwitzHassShortWhiteNoise(params, scene, oi, display);
%%% build outersegment
identityOS = osCreate('identity');
identityOS = osSet(identityOS, 'rgbData', sceneRGB);

%%% build rgc
% rgc1 = rgcLinear(scene, sensor, osIdentity, 'right', 3.75, 180);
% rgc1 = rgcLNP(scene, sensor, osIdentity, 'right', 3.75, 180);
rgc1 = rgcGLM(scene, sensor, identityOS, 'right', 3.0, 180);
% rgc1 = rgcSubunit(scene, sensor, identityOS, 'right', 3.0, 180);

%%% compute rgc
rgc1 = rgcCompute(rgc1, identityOS);

rf1 = (rgc1.mosaic{1}.sRFcenter{1,1} - rgc1.mosaic{1}.sRFsurround{1,1}); rf1 = rf1./max(abs(rf1(:)));
rf2 = (rgc1.mosaic{1}.linearResponse{1,1,2}(:,:,1)); rf2 = rf2./max(abs(rf2(:)));
diff1 = rf1 + rf2;  
max(abs(diff1(:)));

ir1 = (rgc1.mosaic{1}.tCenter{1,1}); ir1 = ir1./max(abs(ir1(:)));
ir2 = (rgc1.mosaic{1}.linearResponse{1,1,1}(:,:,1)); ir2 = ir2./max(abs(ir2(:)));
diff2 = ir1(:) - ir2(:);  
max(abs(diff2(:)));

%%
tolerance = 1e-14;
UnitTest.assertIsZero(max(abs(diff1(:))),'Comparison for spatial RF and spatial impulse response',tolerance);
UnitTest.validationData('spatial RF',rf1);
UnitTest.validationData('spatial impulse response',rf2);

tolerance = 1e-14;
UnitTest.assertIsZero(max(abs(diff2(:))),'Comparison for temporal IR and temporal response to impulse',tolerance);
UnitTest.validationData('temporal IR',rf1);
UnitTest.validationData('temporal response to impulse',rf2);

end

