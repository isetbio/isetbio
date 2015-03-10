function ieManualViewer(mType)
%View ISET manual information located at ImagEval.com
%
%   ieManualViewer(mType)
%
% Open a browser into the ImagEval web-site to show different parts of the
% user manual.
%
% Examples:
%   ieManualViewer
%   ieManualViewer('home')
%   ieManualViewer('application notes')
%   ieManualViewer('iset functions')
%   ieManualViewer('functions')
%
% Copyright ImagEval Consultants, LLC, 2006.

% TODO:  We should permit local reads when the user is off-line.  At
% present, the user has to download the material on his/her own.

% Defaults
if notDefined('mType')
    web('http://imageval.com/documentation/','-browser');
    return;
else   
    mType = ieParamFormat(mType);
    switch mType
        case {'functions','isetfunctions'}
            web('http://www.imageval.com/public/ISET-Functions/ISET','-browser');
        case {'home','imageval'}
            web('http://www.imageval.com','-browser');
        case {'applicationnotes'}
            web('http://imageval.com/application-notes/','-browser');
        otherwise
            web('http://imageval.com/documentation/','-browser');
    end
end

end