function [theRGCmosaicFileNames, prebakedRGCMosaicDir] = listPrebakedMosaics(varargin)
% List the pre-baked mRGC circuits available to this installation
%
% Syntax
%   [theRGCmosaicFileNames, prebakedRGCMosaicDir] = ...
%        mRGCMosaic.listPrebakedMosaics(...)
%
% Optional key/value pairs
%   type - Circuit family to list. Currently only 'ONcenterLinear'.
%
% Returns
%   theRGCmosaicFileNames - Cell array of circuit filenames, in the
%        historical ISETBio naming that loadPrebakedMosaic builds
%   prebakedRGCMosaicDir - The local gallery directory for this family
%
% Description
%   Pre-baked circuits are published in the isetbio-mosaics SDR deposit
%   rather than bundled with the repository. This lists what a user can
%   actually load: any circuit present in the local gallery, plus every
%   circuit available from the deposit. Loading one downloads it on demand.
%
%   Five deposited circuits never had a repository copy, so their legacy
%   filename could not be established by checksum. They are listed with a
%   trailing ' (candidate name, unverified)' marker and cannot be loaded by
%   name until a deposit version publishes their legacy filename.
%
% See also
%   mRGCMosaic.loadPrebakedMosaic, sdrPrebakedCircuitFile,
%   sdrLegacyFilenameMap

    % Parse input
    p = inputParser;
    p.addParameter('type', 'ONcenterLinear', @(x)(ismember(x, {'ONcenterLinear'})));
    p.parse(varargin{:});

    switch (p.Results.type)
        case 'ONcenterLinear'
            rgcTypeSubDirectory = 'ONmRGC';
            collectionName = 'mrgc_on_circuits';
    end

    prebakedRGCMosaicDir = fullfile(isetbioRootPath, 'ganglioncells', 'mosaics', ...
        rgcTypeSubDirectory);

    listing = dir(fullfile(prebakedRGCMosaicDir, '*.mat'));
    localFileNames = {listing.name};

    depositRecords = sdrLegacyFilenameMap(collectionName);

    verifiedNames = {depositRecords.legacy_filename};
    verifiedNames = verifiedNames(~cellfun(@isempty, verifiedNames));

    unresolved = depositRecords(cellfun(@isempty, {depositRecords.legacy_filename}));
    candidateNames = cellfun(@(name) sprintf('%s (candidate name, unverified)', name), ...
        {unresolved.candidate_legacy_filename}, 'UniformOutput', false);

    theRGCmosaicFileNames = unique([localFileNames(:); verifiedNames(:)]);
    theRGCmosaicFileNames = [theRGCmosaicFileNames; candidateNames(:)];
end
