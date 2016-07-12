%
%PAL_PFML_IndependentFit   Determine whether PFs to multiple conditions can
%   be fit individually or whether interdependencies exist.
%
%   syntax: indep = PAL_PFML_IndependentFit(FM)
%
%   Internal function
%
% Introduced: Palamedes version 1.1.0 (NP)
% Modified: Palamedes version 1.4.0 (see History.m)

function indep = PAL_PFML_IndependentFit(FM)

indep = 1;

if isstruct(FM.argsA) || isstruct(FM.argsB) || isstruct(FM.argsG) || isstruct(FM.argsL)
    indep = 0;
end
if ischar(FM.argsA)
    if strncmpi(FM.argsA,'cons',4)
        indep = 0;
    end
end
if ischar(FM.argsB)
    if strncmpi(FM.argsB,'cons',4)
        indep = 0;
    end
end
if ischar(FM.argsG)
    if strncmpi(FM.argsG,'cons',4)
        indep = 0;
    end
end
if ischar(FM.argsL)
    if strncmpi(FM.argsL,'cons',4)
        indep = 0;
    end
end
if isnumeric(FM.argsA)
    if ~isempty(FM.argsA) && ~PAL_isIdentity(FM.argsA)
        indep = 0;
    end
end
if isnumeric(FM.argsB)
    if ~isempty(FM.argsB) && ~PAL_isIdentity(FM.argsB)
        indep = 0;
    end
end
if isnumeric(FM.argsG)
    if ~isempty(FM.argsG) && ~PAL_isIdentity(FM.argsG)
        indep = 0;
    end
end
if isnumeric(FM.argsL)
    if ~isempty(FM.argsL) && ~PAL_isIdentity(FM.argsL)
        indep = 0;
    end
end