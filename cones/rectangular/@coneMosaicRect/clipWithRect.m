function clipWithRect(obj, clippingRect)
% Clip a cone mosaic using a clippingRect
%
% Syntax:
%   clipWithRect(obj, clippingRect)
%
% Description
%    Null all cones outside the clippingRect. The clippingRect must be a
%    structure with fields 'xo', yo', 'width', and 'height', all defined
%    in visual degrees.
%
% Inputs:
%    obj - The cone mosaic object
%    clippingRect - The clipping rect, a struct(see above)
%
% Outputs:
%    None
%
% Optional key/value pairs:
%    None.

% History:
%    2/6/19  NPC  ISETBIO Team, 2019

    % Check that the clippingRect is valid
    coneMosaic.validateClippingRect(clippingRect);
    
    xposDegs = squeeze(obj.patternSupport(:,:,1))*1e6/obj.micronsPerDegree;
    yposDegs = squeeze(obj.patternSupport(:,:,2))*1e6/obj.micronsPerDegree;
    
    idx = xposDegs(:) < clippingRect.xo-clippingRect.width/2;
    obj.pattern(idx) = 0;
    idx = xposDegs(:) > clippingRect.xo+clippingRect.width/2;
    obj.pattern(idx) = 0;
    idx = yposDegs(:) < clippingRect.yo-clippingRect.height/2;
    obj.pattern(idx) = 0;
    idx = yposDegs(:) > clippingRect.yo+clippingRect.height/2;
    obj.pattern(idx) = 0;
end
