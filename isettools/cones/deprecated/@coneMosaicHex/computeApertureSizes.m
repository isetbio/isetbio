function [innerApertureOutlineVarying, outerApertureOutlineVarying, maxApertureMeters, apertureMeters] = ...
    computeApertureSizes(dxInner, dxOuter, innerApertureOutline, ...
    outerApertureOutline, xCoords, yCoords)
% Compute ecc-dependent aperture sizes for use in renderPatchArray
%
% Syntax:
%   [innerApertureOutlineVarying, outerApertureOutlineVarying, maxApertureMeters] = ...
%       computeApertureSizes(...
%       dxInner, dxOuter, innerApertureOutline, outerApertureOutline, ...
%       xCoords, yCoords)
%
% Description:
%     Compute ecc-dependent aperture sizes for use in renderPatchArray
%
% Inputs:
%   dxInner
%   dxOuter, 
%   innerApertureOutline 
%   outerApertureOutline
%   xCoords, 
%   yCoords   
%
% Outputs:
%    innerApertureOutlineVarying - the ecc-varying aperture (inner)
%    outerApertureOutlineVarying - the ecc-varying aperture (outer)
%
% Optional key/value pairs:
%    None.
%        
    xCoords = xCoords(:);
    yCoords = yCoords(:);

    coneEccentricitiesInMeters = (sqrt(xCoords.^2+yCoords.^2));
    coneAnglesInDegrees = atan2(yCoords, xCoords) / pi * 180;

    [~, apertureMeters, ~] = coneSizeReadData(...
        'eccentricity',coneEccentricitiesInMeters,...
        'angle',coneAnglesInDegrees);
    [~,apertureMetersAtZeroEcc, ~] = coneSizeReadData(...
        'eccentricity', 0.0,...
        'angle', 0.0);

    if (~isempty(dxOuter))
        tmp_outerApertureOutline.x = zeros(numel(xCoords), numel(outerApertureOutline.x));
        tmp_outerApertureOutline.y = tmp_outerApertureOutline.x;
    end
    if (~isempty(dxInner))
        tmp_innerApertureOutline.x = zeros(numel(xCoords), numel(innerApertureOutline.x));
        tmp_innerApertureOutline.y = tmp_innerApertureOutline.x;
    end

    increaseFactors = apertureMeters/apertureMetersAtZeroEcc;
    maxApertureMeters = max(increaseFactors)*apertureMetersAtZeroEcc;
    
    for k = 1:numel(xCoords)
        if (~isempty(dxOuter))
            tmp_outerApertureOutline.x(k,:) = outerApertureOutline.x * increaseFactors(k);
            tmp_outerApertureOutline.y(k,:) = outerApertureOutline.y * increaseFactors(k);
        end
        if (~isempty(dxInner))
            tmp_innerApertureOutline.x(k,:) = innerApertureOutline.x * increaseFactors(k);
            tmp_innerApertureOutline.y(k,:) = innerApertureOutline.y * increaseFactors(k);
        end
    end
    if (~isempty(dxOuter))
        outerApertureOutlineVarying.x = tmp_outerApertureOutline.x;
        outerApertureOutlineVarying.y = tmp_outerApertureOutline.y;
    else
        outerApertureOutlineVarying = [];
    end
    if (~isempty(dxInner))
        innerApertureOutlineVarying.x = tmp_innerApertureOutline.x;
        innerApertureOutlineVarying.y = tmp_innerApertureOutline.y;
    else
        innerApertureOutlineVarying = [];
    end
end