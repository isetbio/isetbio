function rootPath = isetbioRootPath()
% Return the path to the root isetbio directory
%
% Syntax:
%   rootPath = isetbioRootPath;
%
% Description:
%    This points at the top level of the isetbio tree on the Matlab path.
%
%    Examples are included within the code.
%
%    This function works by using mfilename to find itself and then walks
%    up from the utility directory to the top level of ISETBio.
%
% Inputs:
%    None.
%
% Outputs:
%    rootPath - The root directory for isetbio
%
% Optional key/value pairs:
%    None.
% 
% See Also:
%    isetbioDataPath, isetRootPath

% History:
%    07/27/17  dhb  Changed to index off of Contents.m, so as not to
%                   clutter isetbio root directory.
%    10/19/17  dhb  Changed to work based on mfilename, the previous method
%                   broke if the root directory was not called 'isetbio'.
%    11/21/17  jnm  Formatting
%    01/16/18  jnm  Formatting update to match Wiki
%    06/12/26  bw   Update walk after promoting utility to the repository root.

% Examples:
%{
    fullfile(isetbioRootPath, 'data')
%}

%% Get path to this function and then walk back up to the isetbio root.
pathToMe = mfilename('fullpath');

%% Walk back up the chain
rootPath = fileparts(fileparts(pathToMe));

%% Older less robust method
% This is an older way of doing it, which breaks when isetbio is not
% installed in a directory called isetbio. This can happen sometimes (if
% you've checked out a branch), so we changed to the above.
% rootPath = which('isetbio/Contents');
% [rootPath,~,~] = fileparts(rootPath);

end
