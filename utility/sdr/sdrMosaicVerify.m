function tf = sdrMosaicVerify(localFile, record)
% True when a local file matches its isetbio-mosaics manifest record
%
% Synopsis
%    tf = sdrMosaicVerify(localFile, record)
%
% Brief description
%   Compare a file's byte count and SHA-256 with the values published in the
%   manifest record. Both must match. This is the single integrity check
%   used for downloaded assets and for cache hits.
%
% Inputs
%   localFile - Full path to the file to check
%   record    - One manifest record from sdrMosaicRecords
%
% Returns
%   tf - Logical scalar
%
% See also
%   sdrMosaicFetch, sdrMosaicRecords, ieHash

details = dir(localFile);
tf = ~isempty(details) && details.bytes == record.bytes && ...
    strcmpi(ieHash(localFile, 'file', 'SHA-256', 'hex'), record.sha256);

end
