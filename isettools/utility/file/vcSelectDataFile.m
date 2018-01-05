function fullName = vcSelectDataFile(dataType, rw, ext, windowTitle)
% Select a data file name for reading or writing
%
% Syntax:
%   fullName = vcSelectDataFile([dataType], [rw], [ext], [windowTitle])
%
% Description:
%    Select a data file name for reading or writing.
%
%    dataType is used to suggest the starting directory. Include stayput
%    (i.e., don't change), or sensor, optics, and any of the directory
%    names inside of data.
%
%    To choose a data file for reading or writing, use this routine. The
%    parameter dataType is a clue about the proper directory to use to find
%    or write the file. 
%
%    To specify whether the file is for reading or writing, use rw = 'r'
%    for reading and rw = 'w' for writing. Default is read.
%
%    You may also pass in an extension to use for filtering file names.
%    Returns fulName = [] on Cancel.
%
% Inputs:
%    dataType    - (Optional) Suggestions of a starting directory
%                  (directory names inside of the data?) Default is
%                  session. Options are:
%                      {'session', 'stayput'}
%                                    - The last selected directory, or the
%                                      current working directory if none
%                                      was selected previously.
%                      {'data'}      - data.
%                      {'algorithm'} - ISET-Algorithms
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
% Notes:
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/29/17  jnm  Formatting

% Examples:
%{
    fullName = vcSelectDataFile()
    fullName = vcSelectDataFile('session', 'r')
    fullName = vcSelectDataFile('session', 'r', 'tif')
    fullName = vcSelectDataFile('data', 'w')
    fullName = vcSelectDataFile('data', 'r')
%}

if notDefined('dataType'), dataType = 'session'; end
if notDefined('rw'), rw = 'r'; end
if notDefined('ext'), ext = '*'; end

curDir = pwd;

% We remember the last directory the user chose. On the first call, this
% variable is empty. But from then on, we use it.
persistent pDir;

switch lower(dataType)
    case {'session', 'stayput'}
        if isempty(pDir)
            fullPath = pwd;
        else
            fullPath = pDir;
        end
    case {'algorithm'}
        fullPath = fullfile(isetRootPath, 'ISET-Algorithms');
        if ~exist(fullPath, 'dir')
            if  ~isempty(pDir)
                fullPath = pDir;
            else
                fullPath = isetRootPath;
            end
        end
    case {'data'}
        fullPath = fullfile(isetbioDataPath);

        % Check that directory exists. If not, try using the last directory
        % Otherwise, just go to data.
        if ~exist(fullPath, 'dir')
            if  ~isempty(pDir)
                fullPath = pDir;
            else
                fullPath = isetRootPath;
            end
        end
    case {'displays'}
        fullPath = fullfile(isetbioDataPath, 'displays');
        if ~exist(fullPath, 'dir')
            if  ~isempty(pDir)
                fullPath = pDir;
            else
                fullPath = isetRootPath;
            end
        end
    otherwise
        if isempty(pDir)
            fullPath = pwd;
        else
            fullPath = pDir;
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