function [idx1, idx2] = ieFieldHeight2Index(fieldHeightList, height)
% Find the field height index closest to a specific height (meters)
%
% Syntax:
%   [idx1, idx2]  = ieFieldHeight2Index(fieldHeightList, height)
%
% Description:
%    You can also request a pair of indices that bound the value. In that
%    case, idx1 < idx2 and 
%
%        fieldHeightList(idx1) <= height <= fieldHeightList(idx2)
%
%    Typically the field heights are in meters and the request is in
%    meters. This routine will run correctly as long as fieldHeightList and
%    height are both in common units. If they are in different units, bad
%    things happen.
%
% Inputs:
%    fieldHeightList - Vector of field heights
%    height          - Height to search for
%
% Outputs:
%    idx1            - Index of fieldHeightList that is equal to height if
%                      only variable queried. Else the index covering the
%                      lower bound if there is no exact match.
%    idx2            - Most often idx1 + 1. If there is no exact match, the
%                      upper bound for height.
%
% Notes:
%    * [Note: XXX - Programming Note: We could return weights that might be
%      used for interpolation]
%    * [Note: JNM - Neither example works - optics not instantiated, at
%      minimum, and I'm not sure what else is missing since
%      rtPSFfieldHeight as a parameter is only used elsewhere in
%      ieParameterOtype, which I'm not entirely positive would tie in
%      correctly here? The only height in opticsGet is 'imageheight']
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    11/30/17  jnm  Formatting

% Examples:
%{
    % NOT WORKING!!!
    % This example is in meters
    optics = opticsCreate;
    fieldHeightList = opticsGet(optics, 'rtPSFfieldHeight');  % meters
    fhIdx = ieFieldHeight2Index(fieldHeightList, 2e-4)
    [idx1, idx2]  = ieFieldHeight2Index(fieldHeightList, 2e-4)
%}
%{
    % NOT WORKING!!!
    % This example is in millimeters
    optics = opticsCreate;
    fieldHeightList = opticsGet(optics, 'rtPSFfieldHeight', 'mm'); 
    fhIdx = ieFieldHeight2Index(fieldHeightList, 0.2)
%}

% This is the index with a value closest to height
[~, idx1] = min(abs(fieldHeightList - height));

% Determine two indices that bound the height value.
if nargout == 2
    if fieldHeightList(idx1) > height
        % Send back the index below. Order everything properly
        idx2 = max(1, idx1 - 1);
        tmp = idx1;
        idx1 = idx2;
        idx2 = tmp;
    else
        % Send back the index above. No need to order
        idx2 = min(length(fieldHeightList), idx1 + 1);
    end
end

end