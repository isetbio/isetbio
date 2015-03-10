function ieInitSession
% Initialize an ISET session file and data
%
%  ieInitSession
%
%  This function initializes the data structures for a basic ISET session.
%  The isetSession file is saved after initialization.
%
% Copyright ImagEval Consultants, LLC, 2005.

global vcSESSION;

% At some point, we will add a name/re-name pull-down to the Session
% window. For now, we simply assign isetSession as the name.
vcSESSION.NAME = sprintf('iset-%s',datestr(now,30));

% Set current directory as the session directory.
vcSESSION.DIR = pwd;

% If we are initiating, then nothing is selected
vcSESSION.SELECTED.SCENE   = [];
vcSESSION.SELECTED.OPTICALIMAGE = [];
vcSESSION.SELECTED.ISA     = [];
vcSESSION.SELECTED.VCIMAGE = [];
vcSESSION.SELECTED.VIDEO   = [];

vcSESSION.SELECTED.GRAPHWIN = [];

% These simulator cell arrays need to be defined to start computing.  For
% now, they are set to empty.  The user sets the parameters through the
% window interface or from the command line.
vcSESSION.SCENE{1}   = [];
vcSESSION.OPTICALIMAGE{1} = [];
vcSESSION.ISA{1}     = [];
vcSESSION.VCIMAGE{1} = [];
vcSESSION.GRAPHWIN   = [];
vcSESSION.VIDEO{1}   = [];

% Start out with the help flag off.
vcSESSION.initHelp = 0;

% Show waitbars or not.
% Make sure the waitbar preference is set.  The others don't slow us down
% much.
iePref = getpref('ISET');
if ~checkfields(iePref,'waitbar')
     setpref('ISET','waitbar',0);
     vcSESSION.GUI.waitbar = 0;
else vcSESSION.GUI.waitbar = iePref.waitbar;
end

% Check gpu computing
try
    gpuArray(0);
    vcSESSION.GPUCOMPUTE = true;
catch 
    vcSESSION.GPUCOMPUTE = false;
end

end