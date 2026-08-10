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
%   ieWebGet, mosaicName, cMosaic

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
    records = localCollectionRecords(collectionName);
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

localFile = localFetchRecord(records(recordIndex));
loadedData = load(localFile, 'cmosaic');
cmosaic = loadedData.cmosaic;

end

function records = localConeMosaicRecords()
records = localCollectionRecords('cone_mosaics');
end

function records = localCollectionRecords(collectionName)
cacheRoot = fullfile(isetbioRootPath, 'data', 'sdr', 'isetbio-mosaics');
topManifest = localFetchManifest(cacheRoot, 'manifest.json');
topData = jsondecode(fileread(topManifest));
collectionIndex = find(strcmp({topData.collections.name}, collectionName), 1);
if isempty(collectionIndex)
    error('mosaicLoad:Manifest', 'The SDR manifest has no %s collection.', collectionName);
end
collectionManifest = localFetchManifest(cacheRoot, topData.collections(collectionIndex).manifest);
collectionData = jsondecode(fileread(collectionManifest));
records = collectionData.files;
end

function localFile = localFetchManifest(cacheRoot, relativePath)
pathParts = split(relativePath, '/');
localFile = fullfile(cacheRoot, pathParts{:});
if ~isfile(localFile)
    ieWebGet('depositname', 'isetbio-mosaics', 'depositfile', relativePath, ...
        'downloaddir', cacheRoot, 'unzip', false);
end
end

function localFile = localFetchRecord(record)
cacheRoot = fullfile(isetbioRootPath, 'data', 'sdr', 'isetbio-mosaics');
pathParts = split(record.path, '/');
localFile = fullfile(cacheRoot, 'cone_mosaics', pathParts{:});
if isfile(localFile) && localRecordMatchesFile(localFile, record)
    return;
end

if isfile(localFile), delete(localFile); end
temporaryRoot = tempname(cacheRoot);
mkdir(temporaryRoot);
cleanup = onCleanup(@() localRemoveFolder(temporaryRoot));
temporaryFile = ieWebGet('depositname', 'isetbio-mosaics', ...
    'depositfile', ['cone_mosaics/' record.path], ...
    'downloaddir', temporaryRoot, 'unzip', false);
if ~localRecordMatchesFile(temporaryFile, record)
    error('mosaicLoad:Checksum', 'Downloaded mosaic failed its manifest checksum or byte-count check.');
end
localDir = fileparts(localFile);
if ~isfolder(localDir), mkdir(localDir); end
movefile(temporaryFile, localFile, 'f');
clear cleanup
end

function tf = localRecordMatchesFile(localFile, record)
details = dir(localFile);
tf = ~isempty(details) && details.bytes == record.bytes && ...
    strcmpi(ieHash(localFile, 'file', 'SHA-256', 'hex'), record.sha256);
end

function localRemoveFolder(folder)
if isfolder(folder), rmdir(folder, 's'); end
end
