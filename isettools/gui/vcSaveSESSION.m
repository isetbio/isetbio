function vcSaveSESSION(fname)
% Save the vcSESSION information
%
% Syntax:
%   vcSaveSESSION(fname)
%
% Description:
%    The vcSESSION structure is saved in the direction vcSESSION.dir. If
%    FNAME is not passed in, or there is no sessionDir, then the user is
%    queried for the file name.
%
% Inputs:
%    fname - String. String containing the session name.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/11/18  jnm  Formatting

global vcSESSION;

if notDefined('fname'), fname = ieSessionGet('name'); end

% Don't save the graph win and GUI figures. These are session dependent.
if checkfields(vcSESSION, 'GRAPHWIN')
    GRAPHWIN = vcSESSION.GRAPHWIN;
else
    GRAPHWIN = [];
end
if checkfields(vcSESSION, 'GUI'), GUI = vcSESSION.GUI; else, GUI = []; end

vcSESSION.GRAPHWIN = [];
vcSESSION.GUI = [];

% Save the session file
sessionDir = ieSessionGet('dir');
if ~exist(sessionDir, 'dir')
    % This can happen if a session directory is changed or placed on a
    % different computer.
    fullSessionFile = vcSelectDataFile('session', 'w');
    sessionDir = fileparts(fullSessionFile);
    ieSessionSet('dir', sessionDir);
else
    fullSessionFile = fullfile(sessionDir, fname);
end

save(fullSessionFile, 'vcSESSION');

% Now, put them back so close can work properly.
vcSESSION.GRAPHWIN = GRAPHWIN;
vcSESSION.GUI = GUI;

fprintf('Saving %s in directory %s\n', fname, sessionDir);

return;
