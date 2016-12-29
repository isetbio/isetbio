function varargout = v_oiTransmittance(varargin)
%
% Validate some of the optical image transmittance calculations

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
ieInit;

%% Human optics
oi = oiCreate;
t = oiGet(oi,'optics transmittance');
w = oiGet(oi,'optics transmittance wave');
if (runTimeParams.generatePlots)
    vcNewGraphWin;
    plot(w,t);
    xlabel('Wave'); ylabel('Transmittance'); grid on
end

%% Create a scene and check that its wavelength is imposed on transmittance
scene = sceneCreate;
scene = sceneSet(scene,'wavelength',450);
oi = oiCompute(oi,scene);
% Notice that the new get respects the new computation
oiGet(oi,'optics transmittance')
if (runTimeParams.generatePlots)
    oiShowImage(oi);
end

%% Diffraction optics
scene = sceneCreate;
scene = sceneSet(scene,'wavelength',550);
oi = oiCreate('diffraction');
oi = oiCompute(oi,scene);
w = 400:10:500;
t = oiGet(oi,'optics transmittance',w);
if (runTimeParams.generatePlots)
    vcNewGraphWin;
    plot(w,t);
    xlabel('Wave'); ylabel('Transmittance'); grid on
    oiShowImage(oi);
end

end
