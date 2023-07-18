function cmosaic = mosaicLoad(sizeDeg,positionDeg)
% Load a stored cMosaic from the data/cones directory
%
% Synopsis
%    cmosaic = mosaicLoad(sizeDeg,positionDeg)
%
% Brief description
%   Load one of the precomputed cone mosaics.  The file format is defined
%   in the routine mosaicName and has the format
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
%   If help, a list of files is returned
%   Otherwise, the cmosaic is returned  (I am not sure if these need
%       updating!)
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
    lst = dir(fullfile(isetbioDataPath,'cones','*.mat'));
    allNames = cell(numel(lst),1);
    fprintf('\n\nSize and Positions available in the cone mosaic library.\n\n')
    fprintf('\n\n            Name  \t\t\t   Size \t   Pos   \n---------------------------------------------------------------------\n')
    for ii=1:numel(lst)
        thisName = lst(ii).name;
        nameParts = split(thisName,'_');
        sz  = cell2mat(textscan(strrep(nameParts{2},'-',' '),'%.1f %.1f')); 
        pos = cell2mat(textscan(strrep(nameParts{3},'-',' '),'%.1f %.1f'));
        fprintf('%s\t\t%02.2f %02.2f\t%02.2f %02.2f\n',thisName,sz,pos);
        allNames{ii} = thisName;
    end
    cmosaic = allNames;
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
        warning('No mosaic exists for that size and position. Here is a list.');
        cmosaic = mosaicLoad('help');
        return;
    end
end

%% Maybe the user sent in both size and position request
fname = mosaicName(sizeDeg,positionDeg);
if exist(fname,'file')
    load(fname,'cmosaic');
else
    warning('No mosaic exists for that size and position. Here is a list.');
    cmosaic = mosaicLoad('help');
end

end