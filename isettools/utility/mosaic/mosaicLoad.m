function cmosaic = mosaicLoad(sizeDeg,positionDeg)
% Load a stored cMosaic from the data/cones directory
%
% Synopsis
%    cmosaic = mosaicLoad(sizeDeg,positionDeg)
%
% Inputs
%   sizeDeg
%   positionDeg
%
% Optional key/val pairs
%
% Outpujt
%
% See also
%   mosaicName, cMosaic

% Examples:
%{
 mosaicLoad('help')
%}
%{
cm = mosaicLoad([1 1],[1 0]);
%}
if isequal(sizeDeg,'help')
    lst = dir(fullfile(isetRootPath,'data','cones','*.mat'));
    fprintf('\n\nSize and Positions available in the cone mosaic library.\n\n')
    fprintf('\n\n  Size   \t   Pos   \n-------------------------\n')
    for ii=1:numel(lst)
        thisName = lst(ii).name;
        nameParts = split(thisName,'_');
        sz  = cell2mat(textscan(strrep(nameParts{2},'-',' '),'%.1f %.1f')); 
        pos = cell2mat(textscan(strrep(nameParts{3},'-',' '),'%.1f %.1f'));
        fprintf('%2.3f %2.3f\t%2.3f %2.3f\n',sz,pos);
    end
    return;
end

%%
fname = mosaicName(sizeDeg,positionDeg);
if exist(fname,'file')
    load(fname,'cmosaic');
else
    error('No mosaic computed for that size and position.')
    % We might offer to compute it and store it now.
end


end