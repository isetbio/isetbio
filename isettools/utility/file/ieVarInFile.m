function bool = ieVarInFile(fullname, varName)
% Check whether a variable is in a particular Matlab file
%
% Syntax:
%   bool = ieVarInFile(fullnameOrVariables, varName)
%
% Description:
%    Check whether a variable is contained in a particular MATLAB file.
%
% Inputs:
%    fullname - Typically, this is the file name. But if you already have
%               the variables loaded via the below
%                     variables = whos('-file', fullname);
%               we notice that fullname is a struct and we treat it as an
%               array of variables.
%    varName  - The variable name string.
%
% Outputs:
%    bool     - The boolean value indicating whether or not the variable is
%               contained in the file.
%
% Notes:
%    * [Note: JNM - Second example does not work the way it appears to be
%      called. Should this be addressed or removed?]

% Examples:
%{
    fullname = fullfile(isetbioDataPath, 'human', 'XYZ.mat');
    ieVarInFile(fullname, 'data')
    ieVarInFile(fullname, 'xyz')
%}
%{
    % [Note: JNM - This function does not perform as what appears to be
    % expected based on the rest of the function. The created variables
    % structure follows 1 row with 9 entries for each of the 4 variables.
    % Calling variables returns the columns from the struct.]

    % Alternate calling convention
    fullname = fullfile(isetbioDataPath, 'human', 'XYZ.mat');
    variables = whos('-file', fullname);
    ieVarInFile(variables, 'xyz')
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