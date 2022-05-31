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

%% Answer a call for help.
if isequal(sizeDeg,'help')
    lst = dir(fullfile(isetRootPath,'data','cones','*.mat'));
    fprintf('\n\nSize and Positions available in the cone mosaic library.\n\n')
    fprintf('\n\n            Name  \t\t\t   Size \t   Pos   \n---------------------------------------------------------------------\n')
    for ii=1:numel(lst)
        thisName = lst(ii).name;
        nameParts = split(thisName,'_');
        sz  = cell2mat(textscan(strrep(nameParts{2},'-',' '),'%.1f %.1f')); 
        pos = cell2mat(textscan(strrep(nameParts{3},'-',' '),'%.1f %.1f'));
        fprintf('%s\t\t%2.3f %2.3f\t%2.3f %2.3f\n',thisName,sz,pos);
    end
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
    else
        error('Could not find %s\n',fname);
    end
end

%% Maybe the user sent in both size and position request
fname = mosaicName(sizeDeg,positionDeg);
if exist(fname,'file')
    load(fname,'cmosaic');
else
    error('No mosaic computed for that size and position.')
    % We might offer to compute it and store it now.
end


end