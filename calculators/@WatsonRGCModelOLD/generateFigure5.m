function hFig = generateFigure5(obj, hFig)
% Generate Figure 5 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigure5();
%
% Description:
%   Generate Figure 5 of Watson (2014) which plots the total RGC RF density
%   as a function of eccentricity for all 4 meridians.
%
% Inputs:
%    obj                       - The WatsonRGCModel object
% Outputs:
%    hFig                      - The generated figure handle
%
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/8/19  NPC, ISETBIO Team     Wrote it.

    figureNumber = '5';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    
    generateRGCRFDensityPlot(obj, @obj.totalRGCRFDensity, obj.eccDegs);
    drawnow;
end