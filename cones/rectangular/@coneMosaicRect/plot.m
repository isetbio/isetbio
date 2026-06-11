function [uData, hf] = plot(obj, plotType, varargin)
% Plot function for @coneMosaicRect class
%
% Syntax:
%   [uData, hf] = plot(obj, plotType, varargin)
%
% Description:
%    There is a specialized plot method for the coneMosaic class that
%    calls this function.
% 
%    When the plot type string begins with 'os ' or 'outersegment ' we pass
%    the arguments along to os.plot(). For example, 
%
%           cMosaic.plot('os impulse response')
%
%    graphs the outer segment impulse response on its own time axis.
%
% Inputs:
%    plotType - string, type of plot (N.B. The '*' means user must select a
%               point on an image to specify the location of the point or
%               line to be plotted.)
%       Cone mosaic              - Color image of the cone arrangement
%
%       Mean absorptions         - Image of the mean absorptions
%       Movie absorptions        - Movie of the absorptions on cone mosaic
%       Mean current             - Image of the mean current
%       Movie current            - Gray scale movie of current
%
%       hline absorptions*       - Graph of horizontal line of absoprtions
%       hline current*           - 
%       hline absorptions lms*   - Three panel graph of LMS absorptions
%       hline current lms*       - Three panel graph of LMS current
%       vline absorptions*       - Vertical line
%       vline current*
%       vline absorptions lms*
%       vline current lms*
%       time series absorptions* - Cone absorptions Graph of a point
%       time series current*     - Cone photocurrent Graph of a point
%
%       Impulse response         - Cone current impulse response
%
%       Cone fundamentals        - Cone pigment without macular or lens
%       Cone spectral QE         - Cone pigment and macular
%       Eye spectral QE          - Cone pigment with macular and lens
%       Macular transmittance    - Graph
%       Macular absorptance      - Graph
%       Macular absorbance       - Graph
%       Eye movement path        - eye movement path
%
%  app - coneRectWindow_App object
%
% Outputs:
%    uData    - The user data
%    hf       - The figure or axes handle
%
% Optional key/value pairs:
%    'hf' - figure handle or control structure, the meaning of value is
%           []: create plot in new figure
%           'none': don't plot
%         (isgraphics figure or axes handle): Use that figure or axes
%    'oi' - Optical image used to get lens transmittance for QE plots 
%    'x' - x value (col) for vline plots 
%    'y' - y value (row) for hline plots
%
% Notes:
%    * TODO: Think about coneImageActivity function at end.
%

% History:
%    xx/xx/16  HJ/BW  ISETBIO TEAM, 2016
%    02/19/18  jnm     Formatting
%
%    I moved all the code that was here into coneRectPlot.  That way
%    we can call that function directly, using coneRectPlot(cm);
%

%  Examples:
%{
  cm = coneMosaicRect;
  cm.plot('cone mosaic');

  % Starts the parallel pool for some reason
  cm.plot('impulse response');
%}

%% Check plot type string if we send this off to the os plot routine

% (Might Find a cleaner way to check and send to os.plot.
%  Maybe create a parse argument string as in ISET.)
if (length(plotType) > 3 && strcmp(plotType(1:3), 'os '))
    obj.os.plot(plotType(4:end), 'cmosaic', obj, varargin{:});
    return;
elseif (length(plotType) > 13 && strcmp(plotType(1:13), 'outersegment '))
    obj.os.plot(plotType(14:end), 'cmosaic', obj, varargin{:});
    return;
end

[uData, hf] = coneRectPlot(obj,plotType,varargin{:});

end


