function bool = ieVarInFile(fullname, varName)
% Check whether a variable is in a particular Matlab file
%
%   bool = ieVarInFile(fullnameOrVariables,varName)
%
% fullname:  Typically, this is the file name.  But if you already have the
%            variables loaded via variables = whos('-file',fullname);, we
%            we notice that fullname is a struct and we treat it as an
%            array of variables.
% varName:   The variable name string.
%
% Example:
%   fullname = fullfile(isetRootPath,'data','human','XYZ.mat');
%   ieVarInFile(fullname,'data')
%   ieVarInFile(fullname,'xyz')
%
%  Alternate calling convention
%
%   variables = whos('-file',fullname);
%   ieVarInFile(variables,varName)
%
% See also:  vcReadImage
%
% Copyright ImagEval Consultants, LLC, 2012.

% Decide if fullname is a file or a list of variables
if ischar(fullname)
    variables = whos('-file', fullname);
elseif isstruct(fullname)
    variables = fullname;
else
    error('unknown fullname format');
end

% Check
bool = ismember(varName, {variables.name});

end