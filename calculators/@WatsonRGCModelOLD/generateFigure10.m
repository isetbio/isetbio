function hFig = generateFigure10(obj, hFig)
% Generate Figure 10 of the Watson 2014 paper
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   WatsonRGCCalc.generateFigures10And11();
%
% Description:
%   Generate Figures 10 of Watson (2014) which plot the midget RGC RF
%   spacing as a function of eccentricity for all 4 meridians.
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

    figureNumber = '10';
    figure(hFig); clf;
    set(hFig, 'Color', obj.figurePrefs.backgroundColor, 'Name', sprintf('Figure %s of %s', figureNumber, obj.paperTitleShort));
    type = 'both polarities';
    generateRGCRFSpacingPlot(obj, @obj.midgetRGCRFSpacing, obj.eccDegs, type);
    title('Midget RGC RF spacing (ON+OFF) as a function of eccentricity');
end