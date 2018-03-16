function oiSum = oiAdd(oi1, oi2, wgts)
% Add the photons in two matched oi structs
%
% Syntax:
%   oiSum = oiAdd(oi1, oi2, wgts)
%
% Description:
%    We use this to form the spectral irradiance images in an @oiSequence. 
%
%       oiSum = oiAdd(oi1, oi2, wgts)
%
%    We should check that the two oi structs match, and then we add the
%    photons from the second oi into the first oi.
%
%    **The function copies everything from oi1, save the photon data
%    information, which is composed of the weighted sums of the two OI
%    structure's photon data.**
%
%    There are examples contained in the code. To access, type 'edit oiAdd'
%    into the Command Window.
%
% Inputs:
%    oi1   - Struct. First OI structure
%    oi2   - Struct. Second OI structure
%    wgts  - Vector. Weight vector for adding the two OI structures.
%
% Outputs:
%    oiSum - Struct. The new OI structure created from the weighted
%            addition of the two inputted OI structures.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - TODO: Fix example. Current has done just enough to kill
%      the error message, there is not an actual illuminance here yet.]
%
% See Also:
%    OISEQUENCE
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    03/01/18  jnm  Formatting

% Examples:
%{
    sc = sceneCreate;
    oi1 = oiCreate;
    illu = oiCalculateIlluminance(oi1);
   oi1 = oiSet(oi1, 'illuminance', illu);
    oi1 = oiAdjustIlluminance(oi1, 10);
    oi2 = oi1;
   wgts = [0.75, 0.25]; 
   oiSum = oiAdd(oi1, oi2, wgts);
%}

oiSum = oi1;
oiSum.data.photons = wgts(1) * oi1.data.photons + ...
    wgts(2) * oi2.data.photons;

end