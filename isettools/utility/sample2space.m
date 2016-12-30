function [r,c] = sample2space(rSamples,cSamples,rowDelta,colDelta)
% Return the physical position of samples
%
%   [r,c] = sample2space(rSamples,cSamples,rowDelta,colDelta)
%
% We treat the center of the samples as (0,0) and use the sampling spacing
% to calculate the locations of the other samples. The locations r and c
% are in the same units as the deltas.
%
% Example 1:
%  [r,c] = sample2space(sceneGet(oi,'rows'),sceneGet(oi,'cols'),sceneGet(oi,'hres'),sceneGet(oi,'wres'));
% 
% Note unfortunate terminology in sceneGet:
%   hres is height resolution; 
%   wres is width resolution
%
% Calculate the positions of every pixel this way:
%    [X,Y] = meshgrid(c,r);
%
% Example 2:
%   cSamples = [1:.2:64];
%   rSamples = [1:.2:64];
%   rowDeltaMicrons = 5;
%   colDeltaMicrons = 5;
%   [rMicrons,cMicrons] = sample2space(rSamples,cSamples,rowDeltaMicrons,colDeltaMicrons)

% Copyright ImagEval Consultants, LLC, 2005.

% Find the center.
%
% Note that although this method 
%   rCenter = max(rSamples(:))/2; 
%   cCenter = max(cSamples(:))/2;
% is faster than what we do, our method correct even if the positions start at 1
rCenter = mean(rSamples(:));
cCenter = mean(cSamples(:));

% Convert to microns relative to center.
r = (rSamples - rCenter)*rowDelta;
c = (cSamples - cCenter)*colDelta;

end

