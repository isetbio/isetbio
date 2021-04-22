% Demo basic usage of the @mRGCMosaic object
%
% Description:
%    Shows basic usage of the new @mRGCmosaic.
%    Here, we generate a midget RGC mosaic object and visualize it along
%    with its input @cMosaic
%
% See Also:
%
%

% History:
%    04/01/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

%% Initialize
ieInit;
clear;
close all;

%% Generate the mosaic
onMRGCmosaic = mRGCMosaic(...
    'whichEye', 'left eye', ...
    'sizeDegs', [.5 .5], ...     % SIZE: 1.0 degs (x) 0.5 degs (y)
    'eccentricityDegs', [0 0] ...  % ECC: (0,0)
    );

xyLims = [-1 1 -1 1];
onMRGCmosaic.visualize(...
    'domainVisualizationLimits', xyLims, ...
    'covisualizeInputConeMosaic', true);

