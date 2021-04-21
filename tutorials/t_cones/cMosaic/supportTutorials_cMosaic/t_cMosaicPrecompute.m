% Demo how to generate a precomputed pair of @cMosaic objects
%
% Description:
%    Shows how to precompute a left/right eye pair of @cMosaic objects. These large
%    mosaics are used to generate smaller mosaics at any eccentricity via
%    cropping.
%

% History:
%    03/23/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

mosaicFOV = 45*[1 1];
maxIterations = 500;

% Generate and save the generation progression of the left eye cone mosaic
cMosaic(...
    'whichEye', 'left eye', ...                     % Generate mosaic for the left eye
    'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', maxIterations, ...         % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', ~true, ...           % visualize the convergence
    'exportMeshConvergenceHistoryToFile', true...
    );

% Generate and save the generation progression of the right eye cone mosaic
cMosaic(...
    'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
    'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', maxIterations, ...         % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', ~true, ...          % visualize the convergence
    'exportMeshConvergenceHistoryToFile', true...
    );

% Extract and export final cone positions for use in all computations
fovDegs = max(mosaicFOV)*1.3;
neuronType = 'cones';
retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'left eye');
retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'right eye');
