classdef rgcLNP < rgcMosaic
% A subclass of rgcMosaic for LNP models of inner retina objects.
%
% Syntax:
%   obj = rgcLNP(rgcL, bipolarM, cellType, [varargin])
%
% Description:
%    rgcMosaic using a LNP (linear-nonlinear-Poisson) computational model
%    rgcLNP is a subclass of rgcMosaic. It is called when creating a new
%    LNP model rgcMosaic for an inner retina object. Typically we get here
%    from the inner retina object with the call:
%
%          ir.mosaicCreate('model', 'LNP', 'type', 'your type goes here')
%
%    The LNP model is detailed in Chichilnisky & Kalmar, J. Neurosci
%    (2002); Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J.
%    Neurosci (2005); and Pillow, Shlens, Paninski, Sher, Litke,
%    Chichilnisky & Simoncelli, Nature (2008).
%
%    The computational model implemented here relies on code by the:
%    <a href="matlab: web(http://pillowlab.princeton.edu/code_GLM.html',...
%    '-browser')">Pillow Lab</a>, which is distributed under the GNU
%    General Public License.
%
%    This function contains examples of usage inline. To access these, type
%    'edit @rgcLNP/rgcLNP.m' into the Command Window.
%
% Inputs:
%    rgcL     - Object. An rgc Layer object.
%    bipolarM - Object. A bipolar mosaic object.
%    cellType - String. A string indicating the cell type.
%
% Outputs:
%    obj      - Object. An rgc mosaic object of model type LNP.
%
% Optional key/value pairs:
%    Needs to be added.
%
% See Also:
%   rgcMosaic.m, rgcLinear.m, rgcGLM.m
%

% History:
%    09/XX/15  JRG  (c) isetbio team
%    07/XX/16  JRG  updated
%    06/20/19  JNM  Documentation pass

% Example:
%{
    ir.mosaicCreate('model', 'LNP', 'type', 'on midget');
%}

%% Properties
    % Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        % [Note: BW - Should this be in rgcMosaic or in each of the special
        % types of rgcMosaics?]

        % NUMBERTRIALS Number of repeats
        numberTrials = 10;

        % Parameter to specify the time bins Pillow uses for coupling and
        % post spike filters (10 ms default)
        dt = 0.1;
    end

    % Protected properties.
    properties (SetAccess = private, GetAccess = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rgcLNP is a subclass of rgcMosaic.
        % See the properties in rgcMosaic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % GENERATORFUNCTION Pillow promotes the linear input voltage using
        % a nonlinear function that he calls the generator function. By
        % default this is an exponential.
        generatorFunction;

        % RESPONSEVOLTAGE The nonlinear voltage response after application
        % of the generator function and the spike coupling responses is
        % represented here
        responseVoltage;

        % POSTSPIKEFILTER These hold the parameters used in the
        % computation. This is the response after a spike
        postSpikeFilter;
    end

    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end

    % Public methods
    methods

        % Constructor
        function obj = rgcLNP(rgcL, bipolarM, cellType, varargin)
            % The constructor for the rgcLNP class
            %
            % Syntax:
            %   obj = rgcLNP(rgcL, bipolarM, cellType, [varargin])
            %
            % Description:
            %    The constructor for the rgcLNP class.
            %
            % Inputs:
            %    rgcL     - Object. An rgc Layer object.
            %    bipolarM - Object. A bipolar mosaic object.
            %    cellType - String. A string indicating the cell type.
            %
            % Outputs:
            %    obj      - Object. a rgc mosaic object of model LNP.
            %
            % Optional key/value pairs:
            %    Needs to be added.
            %

            p = inputParser;
            p.KeepUnmatched = true;

            p.addRequired('rgcL', @(x)(isequal(class(x), 'rgcLayer')));
            p.addRequired('bipolarM', ...
                @(x)(isequal(class(x), 'bipolarMosaic')));
            p.addRequired('cellType', @ischar); % Could check better

            p.parse(rgcL, bipolarM, cellType, varargin{:});

            % Initialize the parent class
            obj = obj@rgcMosaic(rgcL, bipolarM, cellType, varargin{:});

            % Effect of a spike on output voltages
            obj.generatorFunction = @exp;

            % Post spike filter
            obj.postSpikeFilter = ...
                zeros(size(buildPostSpikeFilter(obj.dt)));
        end
    end

    % Methods that must only be implemented (Abstract in parent class).
    methods (Access = public)
        initialize(obj);
    end

    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end

    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
end
