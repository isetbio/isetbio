function t_mRGCMosaicSynthesizeAtStage1(options)
% Denovo synthesis of the spatial position lattices of cones and mRGC RF centers (stage 1)
%
% Syntax:
%   t_mRGCMosaicSynthesizeAtStage1()
%
% Description:
%   Demonstrates how to inspect and/or generate a lattice of the spatial
%   positions of cones (or RFs) in a cone mosaic (or an mRGC RF mosaic)
%
%  This is set up with key/value pairs that demonstate how to select different
%  options. Different choices are illustrated in the examples
%  in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% History:
%    08/28/25  NPC  Wrote it.

% Examples:
%{

% UTTBSkip

% Skip running these examples during autovalidation because some of them load
% a lattice progression file to visualize the generation progress.
% Such files are > 100 MBytes, so they cannot be included in a regular (non-gitlfs)
% github repository.

% NOTE: To run any RGC-related ISETBio code, such as this tutorial, users must follow
% the directions discribed in:
%    https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cell-(RGC)-mosaics
% under section "Configuring ISETBio to access RGC resources and run RGC simulations"
%

    % Example #1: Simply inspect the generation of the default mRGC mosaic lattice
    t_mRGCMosaicSynthesizeAtStage1();

    % Example #2: Synthesize a lattice, here, a lattice of the
    % mRGC mosaic in the right eye which is 16-deg wide
    % Note: synthesizing lattices is an iterative, compute-intense
    % operation which can take many hours to complete
    % depending on the power of the computer it is run on
    t_mRGCMosaicSynthesizeAtStage1(...
        'whichEye', 'right eye', ...
        'neuronType', 'midget ganglion cells', ...
        'sourceLatticeSizeDegs', 16, ...
        'maxIterations', 512, ...
        'onlyInspectLattice', false);

    % Example #3: Inspect the generation of a specific lattice.
    % This example should not be run by regular users as the required lattice
    % progress files exceed the 100 MByte limit of regular github repositories
    % and are therefore not included in ISETBio.
    %
    % This is a multi-step process: Do the following:
    % Step1: Find filenames of lattices that are available with their
    generation progress history included
    theLatticePatchFileNames = retinalattice.listPrecomputedPatches(...
        'withGenerationProgressHistory', true)
    
    % Step2: Pick one of the returned filenames
    theLatticeFileName = theLatticePatchFileNames{1};

    % Step3: Decode the filename to extract the eye, the neuron type, and
    % the source lattice size
    [whichEye, neuronType, sourceLatticeSizeDegs, hasProgress] = ...
        retinalattice.decodeFileName(theLatticeFileName);

    % Step4: Inspect the lattice generation
    t_mRGCMosaicSynthesizeAtStage1(...
        'whichEye', whichEye, ...
        'neuronType', neuronType, ...
        'sourceLatticeSizeDegs', sourceLatticeSizeDegs, ...
        'onlyInspectLattice', true);

%}

arguments
    % --- Which neuron type to generate the spatial position mosaic for ---
    % Currently available types  'cones' and 'midget ganglion cells' 
    options.neuronType (1,:) char {mustBeMember(options.neuronType,{'cones','midget ganglion cells'})} = 'midget ganglion cells';

    % Left or right eye
    options.whichEye (1,:) char {mustBeMember(options.whichEye,{'left eye','right eye'})} = 'right eye';

    % Currently only 64 degs is available
    options.sourceLatticeSizeDegs = 16; 

    % The look up entry size for eccentricity. The higher this is, 
    % the smoother the lattice will be (but will take longer to compute)
    options.eccentricityLookUpTableSamplesNum (1,1) double =  512;

    % How many iterations to terminate the iterative optimization after
    options.maxIterations (1,1) double = 8000;

    % Whether to use parallel processing
    options.useParfor (1,1) logical = true;

    % Whether to export the lattice generation progress
    options.exportHistoryToFile (1,1) logical = true;
    
    % Whether to visualize the lattice position convergence
    options.visualizeConvergence (1,1) logical = true;

    % Whether to only inspect a previously-generated lattice, 
    % not re-generate it
    options.onlyInspectLattice (1,1) logical = true;
end % arguments


% Set the position lattice params
latticeParamsStruct = struct(...
        'neuronType', options.neuronType, ...
        'sourceLatticeSizeDegs', options.sourceLatticeSizeDegs, ...    
        'whichEye', options.whichEye, ...
        'maxIterations', options.maxIterations, ...
        'exportHistoryToFile', options.exportHistoryToFile, ...
        'eccentricityLookUpTableSamplesNum', options.eccentricityLookUpTableSamplesNum, ...
        'visualizeConvergence', options.visualizeConvergence, ... 
        'useParfor', options.useParfor);
 
if (~options.onlyInspectLattice)
    % Generate the lattice
    retinalattice.generatePatch(latticeParamsStruct);
end

% Inspect the lattice
retinalattice.inspectPatch(options.neuronType, options.sourceLatticeSizeDegs, options.whichEye);