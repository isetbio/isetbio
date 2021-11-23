function ofile = save(cmosaic,fname,overwrite)
% Save the cMosaic into a file with the variable name cmosaic
%
% Inputs
%    cmosaic
%    fname  - Can be just a filename or a full path
%    overwrite - What to do about overwriting
%
% Output
%    ofile - Full path to the output file
%
% See also
%

% NOTE:
% We might decide to clear the excitation data before
% saving, if we store the excitation data in here.

if ~exist('overwrite','var'), overwrite = false; end
if ~exist('fname','var'), fname = 'cmosaic.mat'; end

ofile = fname;
[~,~,e] = fileparts(ofile);
if isempty(e), ofile = [ofile,'.mat']; end

if exist(ofile,'file')
    if overwrite
        save(ofile,'cmosaic');
    else
        ofile = uiputfile('*.mat','Select a file',fname);
        save(ofile,'cmosaic');
    end
else
    save(ofile,'cmosaic');
end

% Return Full file path always
ofile = which(ofile);

end