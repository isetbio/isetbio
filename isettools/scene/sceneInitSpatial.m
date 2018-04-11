function scene = sceneInitSpatial(scene)
% Initialize the scene field of view to 10 deg.  
%
% Syntax:
%	scene = sceneInitSpatial(scene)
%
% Description:
%    This field of view is small as most cameras see a 40 deg.  But it is
%    the right size for a small sensor at 8 um and 100x100, as we use for
%    many evaluations.
%
%    We are leaving this trivial routine here for potential future
%    development.
%
%    If the parameter is already set, then it is not modified by this
%    routine.
%
% Inputs:
%    scene - The scene structure
%
% Outputs:
%    scene - The modified scene structure
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/21/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match the Wiki.

% Degrees
if notDefined('scene.wAngular'), scene = sceneSet(scene,'fov', 10); end

end