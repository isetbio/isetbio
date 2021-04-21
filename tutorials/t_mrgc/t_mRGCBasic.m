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
    'whichEye', 'right eye', ...
    'sizeDegs', [1.0 1.5], ...     % SIZE: 1.0 degs (x) 0.5 degs (y)
    'eccentricityDegs', [0 0] ...  % ECC: (0,0)
    );

onMRGCmosaic.inputConeMosaic.visualize()