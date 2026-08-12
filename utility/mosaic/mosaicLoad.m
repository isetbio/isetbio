function cmosaic = mosaicLoad(sizeDeg,positionDeg)
% Load a stored cMosaic from the ISETBio SDR mosaic cache
%
% Synopsis
%    cmosaic = mosaicLoad(sizeDeg,positionDeg)
%    records = mosaicLoad('list')
%    records = mosaicLoad('list',collectionName)
%
% Brief description
%   Load one of the precomputed cone mosaics. Assets are selected from the
%   published isetbio-mosaics manifest, downloaded through ieWebGet, and
%   cached under data/sdr/isetbio-mosaics. The legacy filename format is
%   still accepted.
%
%      sprintf('cmosaic_%.1f-%.1f_%.1f-%.1f.mat',sizeDegs,positionDegs);
%
%       ( cmosaic_rowszdeg-colszdeg_xposdeg-yposdeg )
%
% Inputs
%   sizeDeg
%   positionDeg
%
% Optional key/val pairs
%
% Output
%   If help or list, a manifest-backed list of available files is returned.
%   collectionName may be 'cone_mosaics', 'cone_lattices', 'mrgc_lattices',
%   or 'mrgc_on_circuits'. The default is 'cone_mosaics'.
%   Otherwise, the cmosaic is returned  (I am not sure if these need
%       updating!)
%
% See also
%   ieWebGet, mosaicName, cMosaic, sdrMosaicRecords, sdrMosaicFetch

% Examples:
%{
 mosaicLoad('help')
%}
%{
cm = mosaicLoad([1 1],[1 0]);
%}

%% Answer a call for help or a manifest-backed listing.
if ischar(sizeDeg) && ismember(lower(sizeDeg), {'help', 'list'})
    collectionName = 'cone_mosaics';
    if nargin >= 2 && ischar(positionDeg), collectionName = positionDeg; end
    records = sdrMosaicRecords(collectionName);
    fprintf('\nAvailable SDR assets in %s:\n\n', collectionName);
    for ii = 1:numel(records)
        fprintf('%-5s  %s\n', records(ii).eye, records(ii).path);
    end
    cmosaic = records;
    return;
end

%%  Maybe the user sent in just a name.
if ischar(sizeDeg) 
    % See if it has a .mat extension
    tmp = strsplit(sizeDeg,'.');
    if ~isequal(tmp{end},'mat'), fname = [sizeDeg, '.mat'];
    else, fname = sizeDeg; 
    end
    if exist(fname, 'file')
        load(fname,'cmosaic');
        return;
    end
    bundledFile = fullfile(isetbioDataPath, 'cones', fname);
    if isfile(bundledFile)
        loadedData = load(bundledFile, 'cmosaic');
        cmosaic = loadedData.cmosaic;
        return;
    end

    legacyParts = regexp(fname, '^cmosaic_([0-9.]+)-([0-9.]+)_(-?[0-9.]+)-(-?[0-9.]+)\.mat$', 'tokens', 'once');
    if isempty(legacyParts)
        error('mosaicLoad:UnknownMosaic', ...
            'No local mosaic named ''%s'' and it is not a legacy mosaic filename.', fname);
    end
    sizeDeg = reshape(str2double(legacyParts(1:2)), 1, []);
    positionDeg = reshape(str2double(legacyParts(3:4)), 1, []);
end

%% Maybe the user sent in both size and position request
bundledFile = mosaicName(sizeDeg, positionDeg);
if isfile(bundledFile)
    loadedData = load(bundledFile, 'cmosaic');
    cmosaic = loadedData.cmosaic;
    return;
end

records = localConeMosaicRecords();
recordIndex = find(arrayfun(@(record) ...
    isequal(record.size_degrees(:)', sizeDeg(:)') && ...
    isequal(record.eccentricity_degrees(:)', positionDeg(:)'), records), 1);
if isempty(recordIndex)
    error('mosaicLoad:UnknownMosaic', ...
        'No SDR cone mosaic has size %s at eccentricity %s.', ...
        mat2str(sizeDeg), mat2str(positionDeg));
end

localFile = sdrMosaicFetch('cone_mosaics', records(recordIndex));
loadedData = load(localFile, 'cmosaic');
cmosaic = loadedData.cmosaic;

end

function records = localConeMosaicRecords()
records = sdrMosaicRecords('cone_mosaics');
end
