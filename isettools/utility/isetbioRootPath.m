function rootPath=isetbioRootPath()
%%isetbioRootPath   Return the path to the root isetbio directory
%
% Syntax:
%    rootPath=isetbioRootPath;
%
% Description:
%    This points at the top level of the isetbio tree on the Matlab path.
%    
% Notes:
%    This function works by using the function mfilename to find itself,
%    and then walks back up the result to the top level of isetbio.  Thus,
%    you can't move this function within the isetbio tree without also
%    adjusting the number of levels in the walk to match where you move it to.
% 
% Examples:
%   fullfile(isetbioRootPath,'data')
%
% See also: isetbioDataPath, isetRootPath

% 07/27/17  dhb  Changed to index off of Contents.m, so as not to clutter isetbio root directory.
% 10/19/17  dhb  Changed to work based on mfilename, the previous method broke if the root directory
%                was not called 'isetbio'.

%% Get path to this function and then walk back up to the isetbio root.
pathToMe = mfilename('fullpath');

%% Walk back up the chain
rootPath = fileparts(fileparts(fileparts(pathToMe)));

%% Older less robust method
% This is an older way of doing it, which breaks
% when isetbio is not installed in a directory called
% isetbio.  This can happen sometimes (if you've
% checked out a branch), so we changed to the above.
% rootPath=which('isetbio/Contents');
% [rootPath,~,~]=fileparts(rootPath);

return
