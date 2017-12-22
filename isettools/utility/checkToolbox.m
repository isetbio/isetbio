function [valid, toolboxes] = checkToolbox(toolboxName)
% Checks whether certain matlab toolbox has been installed
%
% Syntax:
%   [valid, toolboxes] = checkToolbox(toolboxName)
%
% Description:
%    Checks whether certain matlab toolbox has been installed
% 
% Inputs:
%    toolboxName - The name of the toolbox you are attempting to verify the
%                  existence of.
%
% Outputs:
%    valid        - The boolean value indicating the installation status of
%                  the desired toolbox.
%    toolboxes    - Array listing all of the toolboxes
%
% Notes:
%

% Examples:
%{
    [valid,tbxList] = checkToolbox('Parallel Computing Toolbox');
%}

% vv contains all of the toolbox names
toolboxes  = ver;
valid = any(strcmp({toolboxes.Name}, toolboxName));

end