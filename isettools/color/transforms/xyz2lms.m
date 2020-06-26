function imgLMS = xyz2lms(imgXYZ, cbType, method, varargin)
% Transform an XYZ image to Stockman LMS cone format, colorblind permitted
%
% Syntax:
%   imgLMS = xyz2lms(imgXYZ, [cbType], [method], [varargin])
%
% Description:
%    This function convert XYZ data to LMS in Stockman space. When cbType
%    is not passed in, or it is set to either 0 or 'trichromats', this is
%    the calculation.
%
%    A calculation for color blind dichromats can also be performed. This
%    is set by using format the cbType variable. In this case, missing cone
%    types are estimated from values of remaining cone types. This is a
%    different from something else one might do for dichromats, namely just
%    return the cone values that are there.
%
%    The estimate for color blind is done either by:
%      * interpolating the missing cone using the algorithm in Brettel,
%        Vienot and Mollon JOSA 14/10 1997. The idea in that paper is that
%        the preserved cones are preserved. The missing cone is assigned a
%        value that is a piece-wise linear transform of the preserved cones
%      * interpolating the missing cone using the linear interpolation
%        proposed by Jiang, Farrell and Wandell, 2015
%      * returning a constant for the missing cone type
%
%    This function contains examples of usage inline. To access these, type
%    'edit xyz2lms.m' into the Command Window.
%
% Inputs:
%    imgXYZ   - Matrix. The XYZ image to transform.
%    cbType   - (Optional) Scalar/String. Type of colorblindness. This can
%               be represented by an integer, or a string, with the
%               following options available:
%      {0, 'trichromat'}             - Default. Trichromatic observer.
%      {1, 'protan', 'protanopia'}   - Protanopia observer, missing L cones
%      {2, 'deutan', 'deuteranopia'} - Deuteranope, missing M cones
%      {3, 'tritan', 'tritanopia'}   - Tritanope, missing S cones
%    method   - (Optional) String. Which algorithm to be used to
%               interpolate the missing cone values. Options (and their
%               varargin information) include:
%        brettel: Default. Use Bettel's algorithm. Provide numeric value to
%                 varargin of a white XYZ value.
%        linear: Use linear interpolation algorithm
%        constant: Leave missing cone values as a constant. Default 0.
%    varargin - (Optional) VARIES. Extra arguments for specific methods.
%               See the method entries above for more information.
%
% Outputs:
%   LMS       - Matrix. The values in Stockman LMS space.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%   * [NOTE - DHB: Following exactly what happens with the varargin values
%      as they are passed down the chain is not trivial. Someday, rewrite
%      with key value pairs so that the code is clearer.]
%   * [NOTE - DHB: Not clear if it is 2 degree or 10 degree XYZ/Stockman.]
%
% See Also:
%   lms2srgb, xyz2lms, colorTransformMatrix
%

% History:
%    xx/xx/15  HJ/BW ISETBIO TEAM Created
%    11/01/17  jnm   Comments, formatting, include code to correctly
%                    re-route 'default' as a cbType entry
%    11/02/17  dhb   I didn't think passing 'default' was a very good thing
%                    to allow. Got rid of it. Might break some calling
%                    code. Added error message that says what to do.
%    11/05/17  dhb   Change 'trichromats' -> 'trichromat'. This is what
%                    the called routine wants.
%    11/17/17  jnm  Formatting
%    07/15/19  JNM  Formatting update

% Examples:
%{
    scene = sceneCreate('reflectance chart');
    vcAddAndSelectObject(scene);
    sceneWindow
    imgXYZ = sceneGet(scene, 'xyz');

    whiteXYZ = sceneGet(scene,'illuminant xyz');
    cbType = 'Tritanope';
    imgLMS = xyz2lms(imgXYZ, cbType, 'Brettel', whiteXYZ);

    vcNewGraphWin;
    imagescRGB(lms2srgb(imgLMS));
%}

%% Check input parameters
if notDefined('imgXYZ'), error('XYZ Image Required'); end
if notDefined('cbType'), cbType = 'trichromat'; end
if notDefined('method'), method = 'Brettel'; end

% if cbType is 'trichomats', reassign to prevent breaking down the line.
if (ischar(cbType) && strcmp(lower(cbType),'default'))
    error('Change passed ''default'' cbType to ''trichromat''');
end
if ischar(cbType) && (strcmp(lower(cbType), 'trichromats') || ...
        strcmp(lower(cbType), 'trichromat'))
    cbType = 0;
end

%% Transform to LMS
%
%  Shape of the inputs
sz = size(imgXYZ);
if numel(imgXYZ) == 3 % a 3-value XYZ vector
    imgXYZ = reshape(imgXYZ, [1 1 3]);
elseif ismatrix(imgXYZ) && size(imgXYZ, 2) == 3
    % Nx3, XW format
    imgXYZ = reshape(imgXYZ, [sz(1) 1 3]);
elseif ndims(imgXYZ) == 3 &&  size(imgXYZ, 3) == 3
    % XYZ Image, do nothing
else
    error('Unknown input XYZ image format');
end

% Transform image to LMS format
xyz2lmsM = colorTransformMatrix('xyz2lms');
imgLMS = imageLinearTransform(imgXYZ, xyz2lmsM);

% Interpolate values of missing cones for colorblind
if strcmp(ieParamFormat(method), 'brettel')
    try
        whiteXYZ = reshape(varargin{1}, [1 1 3]);
        whiteLMS = imageLinearTransform(whiteXYZ, xyz2lmsM);
        imgLMS = lms2lmsDichromat(imgLMS, cbType, method, whiteLMS);
    catch
    end
else
    imgLMS = lms2lmsDichromat(imgLMS, cbType, method, varargin{:});
end

%  Back to input shape
imgLMS = reshape(imgLMS, sz);

end