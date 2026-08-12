function theLatticeFile = sourceLatticeFile(sourceLatticeSizeDegs, neuronType, whichEye)
% Resolve a source lattice to a local file, downloading it when needed
%
% Synopsis
%    theLatticeFile = retinalattice.import.sourceLatticeFile(...
%         sourceLatticeSizeDegs, neuronType, whichEye)
%
% Brief description
%   Source lattices are the retina-wide RF-center position files that
%   cMosaic and mRGCMosaic crop. They are published in the isetbio-mosaics
%   SDR deposit rather than bundled with the repository.
%
%   A lattice generated locally still wins: the repository lattice gallery
%   is searched first, so retinalattice.generatePatch output is used without
%   a network request. Otherwise the matching deposited asset is downloaded
%   into the SDR cache and verified against its manifest checksum.
%
%   Only the final-positions lattices are deposited. The *_progress.mat
%   files written during generation stay in the repository gallery.
%
% Inputs
%   sourceLatticeSizeDegs - Nominal generation field of view, in degrees
%   neuronType            - 'cones' or 'midget ganglion cells'
%   whichEye              - 'left eye' or 'right eye'
%
% Returns
%   theLatticeFile - Full path to a local, verified lattice MAT-file
%
% See also
%   retinalattice.configure, retinalattice.import.finalConePositions,
%   retinalattice.import.finalMRGCPositions, sdrMosaicRecords, sdrMosaicFetch

p = retinalattice.configure(sourceLatticeSizeDegs, neuronType, whichEye);

% A locally generated lattice takes precedence over the deposited copy.
theLatticeFile = fullfile(p.latticeGalleryDir, p.patchFinalPositionsSaveFileName);
if isfile(theLatticeFile), return; end

switch (neuronType)
    case 'cones'
        collectionName = 'cone_lattices';
    case 'midget ganglion cells'
        collectionName = 'mrgc_lattices';
    otherwise
        error('sourceLatticeFile:UnknownNeuronType', ...
            'Unknown neuron type: ''%s''.', neuronType);
end

% The manifest names an eye as 'left' or 'right'.
eyeName = strrep(whichEye, ' eye', '');

records = sdrMosaicRecords(collectionName);
recordIndex = find(arrayfun(@(record) ...
    strcmp(record.eye, eyeName) && ...
    record.generation_fov_degrees == sourceLatticeSizeDegs, records), 1);

if isempty(recordIndex)
    availableList = arrayfun(@(record) ...
        sprintf('%s %d deg', record.eye, record.generation_fov_degrees), ...
        records, 'UniformOutput', false);
    error('sourceLatticeFile:UnknownLattice', ...
        ['No %s source lattice for the %s at %g degrees.\n', ...
        'The deposit provides: %s.'], ...
        strrep(neuronType, ' ', ' '), whichEye, sourceLatticeSizeDegs, ...
        strjoin(availableList, ', '));
end

theLatticeFile = sdrMosaicFetch(collectionName, records(recordIndex));

end
