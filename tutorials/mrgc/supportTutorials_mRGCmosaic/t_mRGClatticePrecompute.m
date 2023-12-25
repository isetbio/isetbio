% Demo how to generate a precomputed pair of @mRGCmosaic objects
%
% Description:
%    Shows how to precompute a left/right eye pair of @mRGCmosaic objects. These large
%    mosaics are used to generate smaller mosaics at any eccentricity via
%    cropping.  More specifically, this stores positions and other related
%    information so that the actual mosaic object at a patch can be
%    computed quickly.  The data are stored in:
%        isetbio/isettools/ganglioncells/data/lattices
%
%    This takes a very long time to run, and I'm not sure we want to write
%    out mosaics by default.  The UTTBSkip comment below prevents this from
%    being autorun.  And I set up a conditional at the bottom so it doesn't
%    write anything out unless you change the default value.

% UTTBSkip

% History:
%    04/21/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

% Configure a conservative parpool manager. This gives at least 8 GB RAM/core
% when the flag is true, and minimizes the chance of a crash because
% of memory fill up.
fancyParpoolControl = true;
if (fancyParpoolControl)
    ASPPManager = AppleSiliconParPoolManager('conservative');
end

% Parameters.
mosaicFOV = 45*[1 1];
maxIterations = 3000;

% Generate and save the generation progression of the left eye mRGC mosaic
mRGCMosaic(...
    'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
    'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', maxIterations, ...         % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', ~true, ...          % visualize the convergence
    'exportMeshConvergenceHistoryToFile', true...
);

% Generate and save the generation progression of the right eye mRGC mosaic
mRGCMosaic(...
    'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
    'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', maxIterations, ...         % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', ~true, ...          % visualize the convergence
    'exportMeshConvergenceHistoryToFile', true...
);

% Extract and export final cone positions for use in all computations.
% The retinalattice.savePositionsAtIteration() method loads the mosaic
% positions at each step during the iterative optimization phase along with 
% a measurement of quality of the lattice at each stop. Then it queries the user
% to select the iteration at which to save the lattice positions. 
% We don't really want to do this write just because someone decided
% to run the tutorial, so set write to false by default.
writeMosaics = false;
if (writeMosaics)
    fovDegs = max(mosaicFOV)*1.3;
    neuronType = 'midget ganglion cells';
    retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'left eye');
    retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'right eye');
end

% Restore previous parpool config
if (fancyParpoolControl)
    ASPPManager.restoreLastParpoolSize();
end
