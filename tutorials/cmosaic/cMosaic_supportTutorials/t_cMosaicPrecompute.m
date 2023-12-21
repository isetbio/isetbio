% Demo how to generate a precomputed pair of @cMosaic objects
%
% Description:
%    Shows how to precompute a left/right eye pair of @cMosaic objects. These large
%    mosaics are used to generate smaller mosaics at any eccentricity via
%    cropping.
%
%    This takes a very long time to run, and I'm not sure we want to write
%    out mosaics by default.  The UTTBSkip comment below prevents this from
%    being autorun.  And I set up a conditional at the bottom so it doesn't
%    write anything out unless you change the default value.

% UTTBSkip

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
    'visualizeMeshConvergence', ~true ...           % visualize the convergence
    );

% Generate and save the generation progression of the right eye cone mosaic
cMosaic(...
    'whichEye', 'right eye', ...                    % Generate mosaic for the right eye
    'sizeDegs', mosaicFOV, ...                      % SIZE: x,y in degs
    'eccentricityDegs', [0 0], ...                  % ECC:  x=0.0 degs, y= 0.0 degs
    'computeMeshFromScratch', true, ...             % generate mesh on-line, will take some time
    'maxMeshIterations', maxIterations, ...         % stop iterative procedure after this many iterations
    'visualizeMeshConvergence', ~true ...          % visualize the convergence
    );

% Extract and export cone positions for use in all computations.
% The retinalattice.savePositionsAtIteration() method loads the mosaic
% positions at each step during the iterative optimization phase along with 
% a measurement of quality of the lattice at each stop. Then it queries the user
% to select the iteration at which to save the lattice positions. 
% We don't really want to do this write just because someone decided
% to run the tutorial, so set write to false by default.
writeMosaics = false;
if (writeMosaics)
    fovDegs = max(mosaicFOV)*1.3;
    neuronType = 'cones';
    retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'left eye');
    retinalattice.savePositionsAtIteration(fovDegs, neuronType, 'right eye');
end