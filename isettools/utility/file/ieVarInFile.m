function bool = ieVarInFile(fullname, varName)
% Check whether a variable is in a particular Matlab file
%
% Syntax:
%   bool = ieVarInFile(fullnameOrVariables, varName)
%
% Description:
%    Check whether a variable is contained in a particular MATLAB file.
%
%    Examples are located within the code. To access the examples, type
%    'edit ieVarInFile.m' into the Command Window.
%
% Inputs:
%    fullname - Typically, this is the file name. But if you already have
%               the variables loaded, say
%                     variables = whos('-file', fullname);
%               we notice the struct and do the right thing.
%
%    varName  - The variable name string.
%
% Outputs:
%    bool     - The boolean value indicating whether or not the variable is
%               contained in the file.
%
% Optional key/value pairs:
%    None.
% 

% Examples:
%{
    fullname = fullfile(isetbioDataPath, 'human', 'XYZ.mat');
    ieVarInFile(fullname, 'data')  % True
    ieVarInFile(fullname, 'XYZ')   % False
%}
%{
    % Alternate calling convention
    fullname = fullfile(isetbioDataPath, 'human', 'XYZ.mat');
    variables = whos('-file', fullname);
    ieVarInFile(variables, 'data')
%}

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