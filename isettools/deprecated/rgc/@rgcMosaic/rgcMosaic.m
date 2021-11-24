classdef rgcMosaic < cellMosaic
% Gnerate an rgcMosaic
%
% Syntax:
%   rgcM = rgcMosaic();
%
% Description:
%    The rgcMosaic class defines a particular RGC type and computational
%    model. The model specifies how to compute the RGC response from the
%    input. The type specifies cell and model parameters.
%
%    The rgcLayer class  holds a collection of rgcMosaics. The mosaics
%    themselves are typically created by a function of the rgcLayer class,
%    such as:
%
%          rgcLayer.mosaicCreate('model', 'LNP', 'type', 'cell type');
%
%    This function contains examples of usage inline. To access these, type
%    'edit rgcMosaic.m' into the Command Window.
%
% Inputs:
%    model - String. The rgcMosaic subclass. Options are: 'LNP', 'GLM'
%    type  - String. The mosaic type. Options are: 'ON Parasol', 'OFF
%            Parasol', 'ON Midget', 'OFF Midget', and 'Small Bistratified'.
%
% Outputs:
%    rgcM  - Object. A rgc mosaic object is created.
%
% Optional key/value pairs:
%    Needs to be added.
%
% Notes:
%    * The RGC models are detailed in the following:
%      Chichilnisky & Kalmar, J. Neurosci (2002); Pillow, Paninski, Uzzell,
%      Simoncelli & Chichilnisky, J. Neurosci (2005); and Pillow, Shlens,
%      Paninski, Sher, Litke, Chichilnisky & Simoncelli, Nature (2008).
%    * The computational model implemented for the coupled GLM model relies
%      on code by <http://pillowlab.princeton.edu/code_GLM.html Pillow>,
%      which is distributed under the GNU General Public License.
%
% References:
%    ISETBIO wiki:
%       <a href="matlab: web(strcat('https://github.com/isetbio/', ...
%       'isetbio/wiki/Retinal-ganglion-cells'), '-browser')">RGCS</a>.
%
% See Also:
%   rgcLNP.m, rgcGLM.m, irCreate, s_initRetina
%

% History:
%    XX/XX/15  JRG/BW  ISETBIO team, 2015
%    06/12/19  JNM     Documentation pass
% Example:
%{
    ir.mosaicCreate('model', 'LNP', 'type', 'on midget');
%}


%% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Was Protected properties.  Changing to Public for debugging, and
% maybe forever
properties (SetAccess = public, GetAccess = public)

    % RFDIAMETER receptive field center diameter in MICRONS
    rfDiameter = [];

    % TCENTER temporal impulse response of the center in dt steps or 1 ms??
    tCenter =[];

    % TSURROUND  and of the surround (1 ms timing by default)
    tSurround = [];

    % TONICDRIVE matrix of baselines drives for linear response; if nonzero,
    % cell spikes occasionally with no input
    tonicDrive;

    % RFDIAMAGNITUDE for making movies of response
    % rfDiaMagnitude;

    % RESPONSELINEAR Store the linear response after convolution
    responseLinear = [];

    % RESPONSESPIKES Store the spike times of the responses
    responseSpikes = [];

    % ELLIPSEMATRIX Store the parameters for the RGC sRF ellipses
    ellipseMatrix = [];
end

properties (Access = public)
end

% Private properties. Only methods of the parent class can set these
properties(Access = private)
end

%% Public methods
methods
    % Constructor
    function obj = rgcMosaic(rgcLayer, bipolarM, cellType, varargin)
        % Initialize an rgcMosaic for a particular cell type
        %
        % Syntax:
        %   rgcMosaic(rgcLayer, cellType, 'input mosaic', val)
        %
        % Description:
        %    rgcMosaic is intialized based on a cell type and the
        %    properties of the bipolar mosaic input. The bipolar mosaic is
        %    specified by its index in the bpLayer object that is attached
        %    to rgcLayer.input.
        %
        % Inputs:
        %    rgcLayer     - Object. A rgcLayer object.
        %    bipolarM     - Struct. A bipolar mosaic structure.
        %    cellType     - String. A string describing the cell type.
        %                   Options are: 'onparasol', 'offparasol',
        %                   'onmidget', 'offmidget', and 'onsbc'.
        %
        % Outputs:
        %    obj          - Object. A rgc mosaic object.
        %
        % Optional key/value pairs:
        %    baselineRate - Numeric. The baseline rate. Default 2.27.
        %

        % History:
        %    XX/XX/17  BW   ISETBIO Team, 2017
        %    06/12/19  JNM  Documentation pass

        p = inputParser;
        p.KeepUnmatched = true;
        p.addRequired('rgcLayer', @(x)(isequal(class(x), 'rgcLayer')));
        p.addRequired('bipolarM', ...
            @(x)(isequal(class(x), 'bipolarMosaic')));
        p.addRequired('cellType', ...
            @(x)(ismember(ieParamFormat(x), obj.validCellTypes)));
        p.addParameter('baselineRate', 2.27, @isscalar);
        p.parse(rgcLayer, bipolarM, cellType, varargin{:});

        % We need the parameters in the layer often enough.
        obj.parent = rgcLayer;

        % So we can build the spatial RFs
        obj.input = bipolarM;

        % This the rgc mosaic type
        obj.cellType = strrep(lower(cellType), ' ', '');

        % Generate spatial RFs of the appropriate size for the cell type
        % and TEE (sets sRFcenter, sRFsurround).
        obj.initSpace(varargin{:});

        % Sets temporal RF properties of tCenter/tSurround
        obj.initTime(rgcLayer);

        %% Initialize the baseline (tonic) drive
        % Tonic drive is the bias or DC term for the linear output of the
        % GLM. If the tonic drive term is greater than 0, then there is a
        % baseline firing rate even when the stimulus input is zero. Units
        % of conditional intensity
        obj.tonicDrive = p.Results.baselineRate * ...
            (ones(size(obj.cellLocation)));
    end

end

properties (Constant)
    % VALIDCELLTYPES
    validCellTypes = {'onparasol', 'offparasol', 'onmidget', ...
        'offmidget', 'onsbc'};
end

%% Methods that must only be implemented (Abstract in parent class).
methods (Access=public)
    function window(obj)
        % Tip: Retrieve guidata using
        %    gui = guidata(obj.figureHandle);
        obj.fig = rgcMosaicWindow(obj);
    end

    function val = timeAxis(obj)
        % Time steps in seconds.  Usually, dt is in 0.1 ms
        val = obj.dt * (1:length(obj.tCenter{1})) * 1e-3;
    end

    % Used to print text in the window
    function str = describe(obj)
        % Describe the RGC mosaic properties
        %
        % Prints the relevant text to a string, which is used in the
        % display window.
        %
        % BW, ISETBIO Team, 2017

        parent = obj.parent;  % Used for size and trials.  Needs help.

        % Cell properties
        str = sprintf('Model: %s\n', class(obj));
        txt = sprintf('Cell type: %s\n', obj.cellType);
        str = addText(str, txt);

        % Mosaic properties
        txt = sprintf('N Trials %d\n', parent.nTrials);
        str = addText(str, txt);
        txt = sprintf('Patch size %d (um)\n', 1e6 * parent.size);
        str = addText(str, txt);

        % Spatial temporal properties
        txt = sprintf('Row, Col: %d, %d\n', size(obj.cellLocation));
        str = addText(str, txt);
        txt = sprintf('Time samples: %d\n', size(obj.responseLinear, 3));
        str = addText(str, txt);

        txt = sprintf('Duration: %.0f ms\n', 1e3 * obj.dt * ...
            size(obj.responseLinear, 3));
        str = addText(str, txt);

    end
end

% Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
end

end
