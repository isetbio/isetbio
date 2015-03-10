function mRGB = macbethPatchData(obj,mLocs,delta,fullData,dataType)
%Return a cell array with the linear RGB values from a vcimage or sensor 
%
%    mRGB = macbethPatchData(obj,mLocs,delta)
%
% Returns the linear RGB values from the sensor or processor window
%
% Example:
%
% See Also:  macbethSelect, vcimageMCCXYZ, macbethColorError
%
% Copyright ImagEval Consultants, LLC, 2011.

if notDefined('obj'),   error('vcimage or sensor required'); end
if notDefined('mLocs'), error('Mid locations required'); end
if notDefined('delta'), error('Patch spacing required'); end
if notDefined('fullData'),fullData = 0; end         % Mean, not all the points
if notDefined('delta'),   dataType = 'result'; end  % Default for vcimage

if fullData  % Every value in the patch
    mRGB = cell(1,24);
    for ii = 1:24
        % mLocs(:,mPatch) is a column vector with (row,col)' for the
        % mPatch.
        theseLocs = macbethROIs(mLocs(:,ii),delta);
        mRGB{ii} = vcGetROIData(obj,theseLocs,dataType);
    end
else  % Mean values from each patch
    mRGB = zeros(24,3);
    for ii = 1:24
        % mLocs(:,mPatch) is a column vector with (row,col)' for the
        % mPatch.
        % This code doesn't work properly for the case of an image sensor.
        % It needs to protect against the NaNs returned in that case.  It
        % works OK for the vcimage.  Fix this some day.
        theseLocs = macbethROIs(mLocs(:,ii),delta);
        mRGB(ii,:) = mean(vcGetROIData(obj,theseLocs,dataType));
    end
end

end