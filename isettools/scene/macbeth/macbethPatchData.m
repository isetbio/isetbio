function mRGB = macbethPatchData(obj, mLocs, delta, fullData, dataType)
% Return a cell array with the linear RGB values from a vcimage or sensor 
%
% Syntax:
%	mRGB = macbethPatchData(obj, mLocs, delta, [fullData], [dataType])
%
% Description:
%    Returns the linear RGB values from the sensor or processor window
%
% Inputs:
%    obj      - The vcimage or sensor to analyze
%    mLocs    - The mid-locations
%    delta    - Patch spacing
%    fullData - (Optional) The mean of the complete data. Default 0.
%    dataType - (Optional) Type of the data to analyze. The default is
%               'results' for a vcimage. Options below are sorted by their
%               image/sensor's type:
%           scene:        {'photons', 'energy', 'reflectance'}
%           opticalimage: {'photons', 'energy'}
%           vcimage:      {'results', 'input'} - Default 'results'
%
% Outputs:
%    mRGB     - The cell array containing the RGB patch data.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Complete the requested fix below (copied here for convenience)
%      "This code doesn't work properly for the case of an image sensor.
%      It needs to protect against the NaNs returned in that case. It works
%      OK for the vcimage. Fix this some day."]
%
% See Also:
%    macbethSelect, vcimageMCCXYZ, macbethColorError
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    01/31/18  jnm  Formatting, change second notDefined for 'delta' to be
%                   for 'dataType' instead. call out TODO formally.

if notDefined('obj'), error('vcimage or sensor required'); end
if notDefined('mLocs'), error('Mid locations required'); end
if notDefined('delta'), error('Patch spacing required'); end
if notDefined('fullData'), fullData = 0; end  % Mean, not all the points
if notDefined('delta'), dataType = 'result'; end  % Default for vcimage

if fullData  % Every value in the patch
    mRGB = cell(1, 24);
    for ii = 1:24
        % mLocs(:, mPatch) is a column vector with the form (row, col)' for
        % the mPatch.
        theseLocs = macbethROIs(mLocs(:, ii), delta);
        mRGB{ii} = vcGetROIData(obj, theseLocs, dataType);
    end
else  % Mean values from each patch
    mRGB = zeros(24, 3);
    for ii = 1:24
        % mLocs(:, mPatch) is a column vector with the form (row, col)' for
        % the mPatch.
        %
        % This code doesn't work properly for the case of an image sensor.
        % It needs to protect against the NaNs returned in that case.  It
        % works OK for the vcimage.  Fix this some day.
        theseLocs = macbethROIs(mLocs(:, ii), delta);
        mRGB(ii, :) = mean(vcGetROIData(obj, theseLocs, dataType));
    end
end

end