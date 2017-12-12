function ieManualViewer(mType)
% View ISET manual information located at ImagEval.com
%
% Syntax:
%   ieManualViewer(mType)
%
% Description:
%    Open a browser into the ImagEval web-site to show different parts of
%    the user manual.
%
% Inputs:
%    mType - Manual Type
%
% Outputs:
%    None.
%
% Notes:
%    * TODO: We should permit local reads when the user is off-line. At
%      present, the user has to download the material on his/her own.
%    * [Note: JNM - The /documentation/ page no longer exists, so I've
%      redirected those calls to point to the imageval home page.]
%

% History:
%    xx/xx/06       Copyright ImagEval Consultants, LLC, 2006.
%    11/22/17  jnm  Formatting & issue addressed in notes section.
%

% Examples:
%{
    ieManualViewer
    ieManualViewer('home')
    ieManualViewer('application notes')
    ieManualViewer('iset functions')
    ieManualViewer('functions')
%}

% Defaults
if notDefined('mType')
    web('http://imageval.com/', '-browser');
    return;
else   
    mType = ieParamFormat(mType);
    switch mType
        case {'functions', 'isetfunctions'}
            web('http://www.imageval.com/public/ISET-Functions/ISET', ...
                '-browser');
        case {'home', 'imageval'}
            web('http://www.imageval.com', '-browser');
        case {'applicationnotes'}
            web('http://imageval.com/application-notes/', '-browser');
        otherwise
            web('http://imageval.com/', '-browser');
    end
end

end