function v_ibio_oiTransmittance(varargin)
%
% Validate some of the optical image transmittance calculations

% varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
%
% end
%
% %% Function implementing the isetbio validation code
% function ValidationFunction(runTimeParams)
% ieInit;

%% Human optics
oi = oiCreate('human');
t = oiGet(oi,'optics transmittance');
w = oiGet(oi,'optics transmittance wave');

ieNewGraphWin;
plot(w,t);
xlabel('Wave'); ylabel('Transmittance'); grid on

%% Create a scene and check that its wavelength is imposed on transmittance
scene = sceneCreate;
wave = 450;
scene = sceneSet(scene,'wavelength',wave);
oi = oiCompute(oi,scene);
oiShowImage(oi);  % Faster than oiWindow()

% Notice that the new get respects the new computation
transmittance = oiGet(oi,'optics transmittance',wave);
assert( numel(transmittance) == numel(oiGet(oi,'wave')) );


%% Diffraction optics
scene = sceneCreate;
wave = 550;
scene = sceneSet(scene,'wavelength',wave);
oi = oiCreate('diffraction');
oi = oiCompute(oi,scene);
oiShowImage(oi);

transmittance = oiGet(oi,'optics transmittance',wave);
assert( numel(transmittance) == numel(oiGet(oi,'wave')) );

%% We still have the transmittance data stored in the optics
wave = 400:10:500;
transmittance = oiGet(oi,'optics transmittance',wave);

ieNewGraphWin;
plot(wave,transmittance,'-o');
xlabel('Wave'); ylabel('Transmittance'); grid on
title('Transmittance')

end
