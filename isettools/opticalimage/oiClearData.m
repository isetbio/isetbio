function oi = oiClearData(oi)
% Clear data from optical image structure
%
% Syntax:
%   oi = oiClearData(oi)
%
% Description:
%    Clear the data from the optical image structure as well as the data in
%    the optics structure within it.
%
%    Clear the optics data is not necessarily a good idea, as these data
%    are mainly the PSF value. The optics are treated as a special case in
%    vcImportObject where the default is to preserve the data. The default
%    is to clear for other objects (sensor, vcimage, and so forth).
%
%    It is also possible to clear the optics data using opticsClearData.
%
% Inputs:
%    oi - Struct. The optical image structure
%
% Outputs:
%    oi - Struct. The modified (cleared) optical image structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * (Copied from below) Do we really want to change the topics here?
%      All this does is clear the OTF data. And there is nothing wrong with
%      that data which is specified by its own wavelength OTF.wave anyway.
%      -- Thoughts?
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/02/18  jnm  Formatting
%    06/25/19  JNM  Minor formatting adjustments

bitDepth = oiGet(oi, 'bit depth');
oi = oiSet(oi, 'data', []);
oi = oiSet(oi, 'bit depth', bitDepth);
oi = oiSet(oi, 'wangular', []);
oi = oiSet(oi, 'depth map', []);

% Do we really want to change the topics here?
% All this does is clear the OTF data. And there is nothing wrong with
% that data which is specified by its own wavelength OTF.wave anyway.

% optics = opticsClearData(oiGet(oi,'optics'));
% oi = oiSet(oi, 'optics', optics);

end