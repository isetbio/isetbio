function rootPath = isetRootPath()
% Return the path to the isetbio isettools directory
%
% Syntax:
%    rootPath = isetRootPath;
%
% Description:
%    Return base bath of the isetbio isettools directory.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    isetbioRootPath, isetbioDataPath.

%% This is easy
rootPath = fullfile(isetbioRootPath, 'isettools', []);

end
