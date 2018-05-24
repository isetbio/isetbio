function mat = ieReadMatrix(defMatrix, fmt, prompt)
% Enter values of a matrix
%
% Syntax:
%   mat = ieReadMatrix(defMatrix, [fmt], [prompt])
%
% Description:
%    The user  types in a set of matrix entries for a matrix of the size of
%    the default matrix, defMatrix. If this is not passed in then, the
%    defMatrix = eye(3).
%
%    The code below contains examples of function usage. To access, type
%    'edit ieReadMatrix.m' into the Command Window.
%
% Inputs:
%    defMatrix - (Optional) Matrix. The default matrix. Default is eye(3),
%                which is [1 0 0; 0 1 0; 0 0 1]
%    fmt       - (Optional) String. Number format. Default is '   %.2e'
%    prompt    - (Optional) String. The string user prompt. Default is
%                'Enter the matrix:'
%
% Outputs:
%    mat       - The generated matrix using the provided format.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/01/18  jnm  Formatting

% Example:
%{
	d = ieReadMatrix(zeros(3, 3));
	d = diag(ieReadMatrix(ones(1, 3), '  %.2f', ' Enter peak [0, 1]'));
%}

if notDefined('defMatrix'), defMatrix = eye(3); end
if notDefined('fmt'), fmt = '   %.2e'; end
if notDefined('prompt'), prompt={'Enter the matrix:'}; end

def = {num2str(defMatrix, fmt)};
dlgTitle = 'Matrix Reader';
lineNo = size(defMatrix, 1);
ReSize = 'on';
answer = inputdlg(prompt, dlgTitle, lineNo, def, ReSize);

if isempty(answer)
    mat = [];
    return;
else
    mat = str2num(answer{1});
end

if (size(mat) ~= size(defMatrix))
    warndlg('Matrix does not match requested size.  Returning null.');
    mat = [];
end

return;
