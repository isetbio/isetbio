function sceneFile = createWorkingFolder(sceneFile, varargin)
% Create a local working folder
%
% Syntax:
%   sceneFile = createWorkingFolder(sceneFile, [varargin])
%
% Description:
%    Because we will be modifying the PBRT file and writing out new
%    resource files for each rendered image, we need a local "working"
%    folder where everything will be stored. This is passed in through
%    'workingDir' or is set by default to isetbiopath/local/
%
% Inputs:
%    sceneFile  - String. The scene file's name at original location
%
% Outputs:
%    sceneFile  - String. The file to write in the working directory
%
% Optional key/value pairs:
%    workingDir - String. The working directory. Default /local inside of
%                 the ISETBio root path.
%

% History:
%    10/13/17  TL   Created
%    12/19/17  jnm  Formatting
%    05/29/19  JNM  Documentation pass

%% Parse inputs
p = inputParser;
p.addRequired('sceneFile', @(x)(exist(sceneFile, 'file')));
addParameter(p, 'workingDir', fullfile(isetbioRootPath, 'local'));
p.parse(sceneFile, varargin{:});
workingDir = p.Results.workingDir;

%% Check to see if workingFolder exists
% If it doesn't exist, let's create it
if(~exist(workingDir, 'dir')), mkdir(workingDir); end

%% Copy scene over to working folder
[path, name, ext] = fileparts(sceneFile);
if isempty(path), error('Scene file must be an absolute path.'); end
status = copyfile(path, workingDir);
if ~status
    error('Could not copy scene directory to working directory.');
else
    fprintf('Copied directory from:\n');
    fprintf('%s \n', path);
    fprintf('to \n');
    fprintf('%s \n \n', workingDir);
end

sceneFile = fullfile(workingDir, strcat(name, ext));

end