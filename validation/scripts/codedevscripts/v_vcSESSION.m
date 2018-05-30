function varargout = v_vcSESSION(varargin)
%
% Tests the ieAddObject function
%
% Copyright Imageval LLC, 2018

%{
v_vcSESSION;
%}
varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Initialize ISETBIO
ieInit;
scene = sceneCreate;
ieAddObject(scene);
tst = vcSESSION.SCENE{1};

UnitTest.assert(isequal(tst,vcGetObject('scene')),'ieAddObject/vcGetObject pass for scene');

%%
oi = oiCreate;
oi = oiCompute(oi,scene);

ieAddObject(oi);
tst = vcSESSION.OPTICALIMAGE{1};
UnitTest.assert(isequal(tst,vcGetObject('oi')),'ieAddObject/vcGetObject pass for oi');

end
