function fullName = vcSelectDataFile(dataType, rw, ext, windowTitle)
% Select a data file name for reading or writing
%
% Syntax:
%   fullName = vcSelectDataFile([dataType], [rw], [ext], [windowTitle])
%
% Description:
%    Select a data file name for reading or writing.  This routine
%    maintains a persistent variable that remembers your last choice.
%    This makes it simpler for the user to open the selection in the
%    same directory multiple times, without having to change the
%    command window's directory.
%
%    dataType is used to suggest the starting directory. There are
%    various shortcut options (see below), or you could use a full
%    directory path.
%
%    To specify whether the file is for reading or writing, use rw = 'r'
%    for reading and rw = 'w' for writing. Default is read.
%
%    You may also pass in an extension to use for filtering file names.
%    Returns fullName = [] on Cancel.
%
%    Examples are located within the code. To access the examples, type
%    'edit vcSelectDataFile.m' into the Command Window.
%
% Inputs:
%    dataType    - (Optional) Suggestions of a starting directory. Default
%                  is ''. Options are:
%                      {'', 'stayput'} - Default. The last selected
%                                        directory, stored in the
%                                        persistent directory variable.
%                      {'data'}        - isetbioDataPath
%                      {'sub-directory of data'}
%                                      - A string that is one of the
%                                        sub-directories of isetbioDataPath
%                                        'bipolar', 'color', 'fonts',
%                                        'lights', 'images', 'optics',
%                                        'pbrtscenes', 'surfaces'
%                      {directory-string}
%                                      - Name of an existing directory
%    rw          - (Optional) Read/Write determination of the file. Default
%                  is 'r'. Options are 'r' for read, and 'w' for write.
%    ext         - (Optional) A file extension used to filter filenames.
%                  Default is '*'.
%    windowTitle - (Optional) String for the read window to help the user
%                  know the desired purpose. Default depends on rw. If
%                  read, the default is 'ISET: Read Data', if write, the
%                  default is 'ISET: Write Data'
%
% Outputs:
%    fullName    - The full file and path name of the data file. If the
%                  operation is canceled, will be [].
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/29/17  jnm  Formatting
%    01/29/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip.  Requires user input.
    fullName = vcSelectDataFile()
    fullName = vcSelectDataFile('', 'r')
    fullName = vcSelectDataFile('', 'r', 'tif')
    fullName = vcSelectDataFile('data', 'w')
    fullName = vcSelectDataFile('color', 'r')
%}

%%
if notDefined('dataType'), dataType = ''; end
if notDefined('rw'), rw = 'r'; end
if notDefined('ext'), ext = '*'; end

curDir = pwd;

% We remember the last directory the user chose. On the first call, this
% variable is empty. But from then on, we use it.
persistent pDir;

%%
switch lower(dataType)
    case {'stayput', ''}
        % Use the persistent directory name we have stored
        if isempty(pDir),   fullPath = pwd;
        else,               fullPath = pDir;
        end
    case {'data'}
        % Go to the isetbio data directory
        fullPath = fullfile(isetbioDataPath);

    case {'bipolar', 'color', 'fonts', 'lights', 'images', 'optics', ...
            'pbrtscenes', 'surfaces'}
        fullPath = fullfile(isetbioDataPath, dataType);
  
    otherwise
        if exist(dataType, 'dir'), fullPath = dataType;
        else
            error('Could not find directory %s\n', dataType);
        end
end

chdir(fullPath);
fileFilter = ['*.', ext];
switch lower(rw)
    case 'r'
        if notDefined('windowTitle'), windowTitle = 'ISET: Read Data'; end
        [fname, pname] = uigetfile(fileFilter, windowTitle);
    case 'w'
        if notDefined('windowTitle'), windowTitle = 'ISET: Write Data'; end
        [fname, pname] = uiputfile(fileFilter, windowTitle);
    otherwise
        error('Read/Write set incorrectly')
end

% Clean up and return
chdir(curDir)
if isequal(fname, 0) || isequal(pname, 0)
    fullName = [];
    disp('User canceled');
else
    fullName = fullfile(pname, fname);
    pDir = pname;
end

end