function imgLMS = xyz2lms(imgXYZ, cbType, method, varargin)
% Transform an XYZ image to Stockman cone format, colorblind permitted
%
% Syntax:
%   imgLMS = xyz2lms(imgXYZ, cbType, varargin)
%
% Description:
%    This function convert XYZ data to LMS in Stockman space. When cbType
%    is not passed in, or it is set to either 0, 'default' or
%    'trichromats', this is the default calculation.
%
%    A calculation for color blind can also be performed.  This is set by
%    using format the cbType variable. The estimate for color blind is done
%    either by:
%      * interpolating the missing cone using the algorithm in Brettel,
%        Vienot and Mollon JOSA 14/10 1997. The idea in that paper is that
%        the preserved cones are preserved. The missing cone is assigned a
%        value that is a piece-wise linear transform of the preserved cones
%      * interpolation the missing cone using the linear interpolation
%        proposed by Jiang, Joyce and Wandell, 2015
%      * returning a zero for the missing cone type  
%
% Inputs:
%    imgXYZ   - XYZ image to transform
%    cbType   - Type of colorblindness, can be choosen from
%      {0, 'Trichromats'}            - Trichromatic observer (Default)
%      {1, 'Protan', 'Protanopia'}   - Protanopia observer, missing L cones
%      {2, 'Deutan', 'Deuteranopia'} - Deuteranope, missing M cones
%      {3, 'Tritan', 'Tritanopia'}   - Tritanope, missing S cones
%
%    method   - Algorithm to be used to interpolate the missing cone values
%      'Brettel'  - Using Bettel's algorithm (Default)
%      'Linear'   - Using Linear Interpolation
%      'Constant' - leave missing cone values as a constant (default 0)
%
% Optional Key/Value Pairs:
%      varargin{1} - white XYZ values for Brettel method
%                    extrapolate values for Constant method
%
% Outputs:
%   LMS       - values in Stockman LMS space
%
% See also: 
%   lms2srgb, xyz2lms, colorTransformMatrix
%
% HJ/BW, ISETBIO TEAM, 2015

% Examples:
%{
   scene = sceneCreate('reflectance chart');
   vcAddAndSelectObject(scene); sceneWindow
   imgXYZ = sceneGet(scene, 'xyz');
   
   whiteXYZ = sceneGet(scene,'illuminant xyz');
   cbType = 'Tritanope'; 
   imgLMS = xyz2lms(imgXYZ, cbType, 'Brettel', whiteXYZ);
   
   vcNewGraphWin; imagescRGB(lms2srgb(imgLMS));
%}

%% Check input parameters
if notDefined('imgXYZ'), error('XYZ Image Required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('method'), method = 'Brettel'; end

%% Transform to LMS
%  Shape of the inputs
sz = size(imgXYZ);
if numel(imgXYZ) == 3 % a 3-value XYZ vector
    imgXYZ = reshape(imgXYZ, [1 1 3]);
elseif ismatrix(imgXYZ) && size(imgXYZ, 2) == 3 % Nx3, XW format
    imgXYZ = reshape(imgXYZ, [sz(1) 1 3]);
elseif ndims(imgXYZ) == 3 &&  size(imgXYZ, 3) == 3 % XYZ Image, do nothing
else error('Unknown input XYZ image format'); 
end

% transform image to LMS format
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