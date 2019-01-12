% Show how to add and get objects from the vcSESSION data base
%
% Syntax:
%   t_codeSESSION
%
% Description:
%    Tests whether the ieAddObject and vcGetObject routines are running
%    properly, interacting with the vcSESSION global variable.
%
%    The vcSESSION global variable registers different scenes and
%    optical images.  The scene and oi structs are stored in vcSESSION.
%    Usually the values of global variables are set or read using
%    ieSessionGet or ieSessionSet.
%
% See Also:
%    ieSessionSet, ieSessionGet, ieAddObject, vcGetObject
%

%% Initialize the database variable
ieInit

%%
scene = sceneCreate;
ieAddObject(scene);
tst = vcSESSION.SCENE{1};

if ~isequal(tst, vcGetObject('scene'))
    warning('scene add object failed');
end

%%
oi = oiCreate;
oi = oiCompute(oi, scene);

ieAddObject(oi);
tst = vcSESSION.OPTICALIMAGE{1};

if ~isequal(tst, vcGetObject('oi'))
    warning('scene add object failed');
end

%%