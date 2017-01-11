function FigureSave(figName,figHandle,figType)
% FigureSave(figName,figHandle,figType)
%
% Encapsulates figure saving.  This way you
% can swap in the method here.
%
% Loops over list in figType if figType is a cell array.
%
% If version is 8 (2014b) or later, forces 'eps' -> 'epsc'
% so that color is preserved in eps output.
%
% 4/6/14  dhb  Wrote it
% 7/14/14 dhb  savefig -> savefigghost
% 8/1/14  dhb  Use saveas for 2014b and after, because of new figure handle
%              stuff.
% 8/3/14  dhb  Fix eps color issue in 2014b by forcing type -> epsc.
%              Trying to get Matlab native to crop output properly.
% 11/14/14 npc Removed {i} from figType{i} in line 53 which was causing a crash.
% 1/6/14  dhb  Trying to get something to work under 2014b, added export_fig

% Set figure itself to be white. This might be annoying for someone who
% doesn't want it white, but we are trying to reproduce the pre-2014b
% default behavior.
set(figHandle,'Color','white');

if iscell(figType)
    for i = 1:length(figType)
        if (exist('export_fig','file'))
            export_fig(figName,['-' figType{i}],figHandle,['-p20']);
        elseif (verLessThan('matlab','8.4'))
            if (exist('exportfig','file'))
                exportfig(figHandle,figName,'FontMode', 'fixed','FontSize', 12,'Width',6,'Height',6, 'color', 'cmyk','Format',figType{i});
            elseif (ismac & exist('savefigghost','file'))
                savefigghost(figName,figHandle,figType{i});
            else
                saveas(figHandle,figName,figType{i});
            end
        else
            if (strcmp(figType{i},'eps'))
                figType{i} = 'epsc';
            end
            if (strcmp(figType{i},'epsc') || strcmp(figType{i},'pdf'))
                print(figHandle,['-d' figType{i}],'-loose',figName);
            else
                saveas(figHandle,figName,figType{i});
            end
        end
    end
else
    if (exist('export_fig','file'))
        export_fig(figName,['-' figType],figHandle,['-p20']);
    elseif (verLessThan('matlab','8.4'))
            if (exist('exportfig','file'))
                exportfig(figHandle,figName,'FontMode', 'fixed','FontSize', 12,'Width',6,'Height',6, 'color', 'cmyk','Format',figType);
            elseif (ismac & exist('savefigghost','file'))
                savefigghost(figName,figHandle,figType);
            else
                saveas(figHandle,figName,figType);
            end
        else
            if (strcmp(figType,'eps'))
                figType = 'epsc';
            end
            if (strcmp(figType,'epsc') || strcmp(figType,'pdf'))
                print(figHandle,['-d' figType],'-loose',figName);
            else
                saveas(figHandle,figName,figType);
            end
        end
    end
