function demoWatsonRGCModel
% Demo different ways of using the WatsonRGCModel class
%
% Syntax:
%   demoWatsonRGCModel();
%
% Description:
%  Run the model to reproduce various figures of the Watson 2014 paper.
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    obj = WatsonRGCModel();
    plotlabOBJ = obj.setUpPlotLab();
    
    WatsonRGCModel.unitTestFigure1('plotlabOBJ', plotlabOBJ);
    WatsonRGCModel.unitTestFigure5('plotlabOBJ', plotlabOBJ);
    WatsonRGCModel.unitTestFigure9('plotlabOBJ', plotlabOBJ);
    WatsonRGCModel.unitTestFigure10('plotlabOBJ', plotlabOBJ);
    WatsonRGCModel.unitTestFigure11('plotlabOBJ', plotlabOBJ);
    WatsonRGCModel.unitTestFigure14('plotlabOBJ', plotlabOBJ);
end

