function t_mRGCMosaicSynthesizeAtStage1(options)
% Denovo synthesis of the spatial position lattices of cones and mRGC RF centers (stage 1)
%
% Syntax:
%   t_mRGCMosaicSynthesizeAtStage1();
%
% Description:
%   Demonstrates how to generate an spatial position lattices of the input
%   cone mosaic and/or the mRGC RF mosaic
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
    t_mRGCMosaicSynthesizeAtStage1();
%}

arguments
    % --- Which neuron type to generate the spatial position mosaic for ---
    % Currently available types  'cones' and 'midget ganglion cells' 
    options.neuronType (1,:) char {mustBeMember(options.neuronType,{'cones','midget ganglion cells'})} = 'midget ganglion cells';

    % Left or right eye
    options.whichEye (1,:) char {mustBeMember(options.whichEye,{'left eye','right eye'})} = 'right eye';

    % Currently only 64 degs is available
    options.sourceLatticeSizeDegs = 64; 

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