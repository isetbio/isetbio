function [theLatticePatchFileNames, latticeGalleryDir] = listPrecomputedPatches(varargin)
% List the midget-RGC source lattices available to this installation
%
% Syntax
%   [theLatticePatchFileNames, latticeGalleryDir] = ...
%        retinalattice.listPrecomputedPatches(...)
%
% Optional key/value pairs
%   withGenerationProgressHistory - When true, list the lattice generation
%        progress files instead of the final-position lattices. Progress
%        files are large generation artifacts that are neither shipped with
%        ISETBio nor published in the SDR deposit, so this list is empty
%        unless lattices were generated locally.
%
% Returns
%   theLatticePatchFileNames - Cell array of lattice filenames, in the
%        historical ISETBio naming that retinalattice.decodeFileName parses
%   latticeGalleryDir - The local lattice gallery directory
%
% Description
%   Final-position lattices are published in the isetbio-mosaics SDR deposit
%   rather than bundled with the repository. This function lists what a user
%   can actually load: any lattice present in the local gallery, plus every
%   lattice available from the deposit. A locally generated lattice is listed
%   once, even when the deposit publishes the same name.
%
% See also
%   retinalattice.decodeFileName, retinalattice.import.sourceLatticeFile,
%   sdrLegacyFilenameMap

    % Parse input
    p = inputParser;
    p.addParameter('withGenerationProgressHistory', false, @islogical);
    p.parse(varargin{:});

    withGenerationProgressHistory = p.Results.withGenerationProgressHistory;

    latticeGalleryDir = fullfile(isetbioDataPath, 'rgc', 'lattices');

    listing = dir(fullfile(latticeGalleryDir, '*.mat'));
    localFileNames = {listing.name};

    if (withGenerationProgressHistory)
        % Generation artifacts are local only.
        theLatticePatchFileNames = localFileNames(...
            contains(localFileNames, 'progress'));
    else
        localFinalFileNames = localFileNames(...
            ~contains(localFileNames, 'progress'));

        % Add the lattices published in the deposit, named as ISETBio has
        % always named them so decodeFileName still parses the result.
        depositRecords = sdrLegacyFilenameMap('mrgc_lattices');
        depositFileNames = {depositRecords.legacy_filename};
        depositFileNames = depositFileNames(~cellfun(@isempty, depositFileNames));

        theLatticePatchFileNames = unique(...
            [localFinalFileNames(:); depositFileNames(:)]);
    end

    theLatticePatchFileNames = theLatticePatchFileNames(:);
end
