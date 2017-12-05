function ret = checkToolbox(toolboxName)
% Checks whether certain matlab toolbox has been installed
%
% Syntax:
%   ret = checkToolbox(toolboxName)
%
% Description:
%    Checks whether certain matlab toolbox has been installed
% 
% Inputs:
%    toolboxName - The name of the toolbox you are attempting to verify the
%                  existence of.
%
% Outputs:
%    ret         - The boolean value indicating the installation status of
%                  the desired toolbox.
%
% Notes:
%    * [Note: JNM - The function was marked incorrectly below as returning
%      all of the toolbox names. It is possible to do that if you change
%      ret above to [ret, vv], but as it stands, the function currently
%      returns a boolean indicating whether or not the specific toolbox is
%      installed on the local machine.]
%

% Examples:
%{
    ret = checkToolbox('Parallel Computing Toolbox');
%}

% vv contains all of the toolbox names
vv  = ver;
ret = any(strcmp({vv.Name}, toolboxName));

end