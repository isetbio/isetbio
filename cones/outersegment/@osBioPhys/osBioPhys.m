classdef osBioPhys < outerSegment 
% Create a biophysically based outersegment (os) object.
%
% Syntax:
%   os = osBioPhys;
%   os = osBioPhys('eccentricity', 0);
%
% Description:
%    This class provides methods and parameters for converting a movie of
%    isomerizations (R*) to outer segment current (pA).
% 
%    Rieke and colleagues defined a set of difference equations as a
%    simulation of the phototransduction cascade. This object defines the
%    parameters and methods to transform the computed isomerizations (R*)
%    in the coneMosaic current.
%
%    If the noiseFlag property of the osLinear object is true, this method
%    adds noise to the current output signal.
%
%    The osBioPhys model is also the basis of how we find the linear
%    filters in the osLinear model, another subclass of outerSegment.
%
%    At present, we just have two sets of parameters, foveal and
%    peripheral. One day we might try to handle eccentricity more finely.
%
% Inputs:
%    None required.
%
% Outputs:
%    The created outersegment object.
%
% Optional key/value pairs:
%    'eccentricity' - Eccentricity in degrees. Determines parameters used.
%                     Currently we just have foveal and peripheral
%                     parameters, and somewhat arbitrarily set the cuttoff
%                     at 10 degrees. Default is 15 degrees.
%
% References:
%    http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%    https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%

% History:
%    xx/xx/16  JRG/HJ/BW  ISETBIO Team, 2016
%
%    08/05/17   dhb         Add eccentricityDegs key/value pair, and comments
%                           about it. This had been suggested previously in
%                           the comments as required. I did it in a backwards
%                           compatiable fashion.
%               dhb         Deprecate 'osType' and change 'eccentricityDegs'
%                           to 'eccentricity'.
%    02/14/18   jnm         Formatting
%    10/19/2020 npc         Made fovealPeripheralCutoffDegs a constant

    properties (Constant)
         %FOVEALPERIPHERALCUTOFFDEGS  Eccentricity in degrees beyond which
        %   we switch from foveal to peripheral parameters.
        %
        %   The 10 degrees value comes from Fred Rieke. Soon we will have
        %   intermediate eccentricity dynamics (3-5 degrees) so this will
        %   change.
        fovealPeripheralCutoffDegs = 10;
    end
    
    properties(Access = private)
        %STATE  Biophysics parameter state
        state; 
    end

    properties(SetAccess = protected, GetAccess = public)
        %MODEL  Structure with biophysical model parameters
        model;
    end

    methods
        % Constructor
        function obj = osBioPhys(varargin)
            % Initialize the osBioPhys
            %
            % Syntax:
            %   obj = osBioPhys([varargin])
            %
            % Description:
            %    Initialize the osBioPhys
            %
            % Inputs:
            %    None required.
            %
            % Outputs:
            %    obj - The created outersegment bioPhysical object
            %
            % Optional key/value pairs:
            %    None.
            %

            % Parse input
            p = inputParser;
            addParameter(p, 'eccentricity', 15, @isnumeric);
            p.parse(varargin{:});
            obj.eccentricityDegs = p.Results.eccentricity; 

            % If eccentricity is greater than foveal/peripheral cutoff, use
            % peripheral parameters.
            if (obj.eccentricityDegs > obj.fovealPeripheralCutoffDegs)
                    %% Peripheral parameters
                    %
                    % sigma     - rhodopsin activity decay rate (1/s)
                    % phi       - phosphodiesterase activity decay rate 1/s
                    % eta       - phosphodiesterase activation rate
                    %             constant (1/s).
                    % gdark     - concentration of cGMP in darkness
                    % k         - constant relating cGMP to current
                    % h         - cooperativity for cGMP->current
                    % cdark     - dark calcium concentration
                    % beta      - rate constant for calcium removal in 1/s
                    % betaSlow  - rate constant for slow calcium modulation
                    %             of channels.
                    % n         - cooperativity for cyclase, hill coef
                    % kGc       - hill affinity for cyclase
                    % OpsinGain - So stimulus can be in R*/s (rate of
                    %             increase in opsin activity per R*/s)
                    obj.model.sigma = 22;       % Default 22
                    obj.model.phi = 22;         % Default 22
                    obj.model.eta = 2000;       % Default 2000
                    obj.model.gdark = 20.5;     % Default 20.5
                    obj.model.k = 0.02;         % Default 0.02
                    obj.model.h = 3;            % Default 3
                    obj.model.cdark = 1;        % Default 1
                    obj.model.beta = 9;           % Default 9
                    obj.model.betaSlow = 0.4;   % Default 0.4
                    obj.model.n = 4;             % Default 4
                    obj.model.kGc = 0.5;        % Default 0.5
                    obj.model.OpsinGain = 10;   % Default 10

            % Othewise use foveal parameters
            else
                    %% Foveal parameters
                    %
                    % sigma     - rhodopsin activity decay rate (1/s)
                    % phi       - phosphodiesterase activity decay rate 1/s
                    % eta       - phosphodiesterase activation rate
                    %             constant (1/s).
                    % gdark     - concentration of cGMP in darkness
                    % k         - constant relating cGMP to current
                    % h         - cooperativity for cGMP->current
                    % cdark     - dark calcium concentration
                    % beta      - rate constant for calcium removal in 1/s
                    % betaSlow  - rate constant for slow calcium modulation
                    %             of channels.
                    % n         - cooperativity for cyclase, hill coef
                    % kGc       - hill affinity for cyclase
                    % OpsinGain - So stimulus can be in R*/s (rate of
                    %             increase in opsin activity per R*/s)
                    obj.model.sigma = 10;       % Default 22
                    obj.model.phi = 22;         % Default 22
                    obj.model.eta = 700;        % Default 2000
                    obj.model.gdark = 20.5;     % Default 20.5
                    obj.model.k = 0.02;         % Default 0.02
                    obj.model.h = 3;            % Default 3
                    obj.model.cdark = 1;        % Default 1
                    obj.model.beta = 5;         % Default 9
                    obj.model.betaSlow = 0.4;   % Default 0.4
                    obj.model.n = 4;            % Default 4
                    obj.model.kGc = 0.5;        % Default 0.5
                    obj.model.OpsinGain = 12;   % Default 10
            end

        end

        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end

        function val = get(obj, varargin)
            val = osGet(obj, varargin{:});
        end

        state = osAdaptSteadyState(obj, bgR, varargin);

        [adaptedData, state] = osAdaptTemporal(pRate, obj);

    end

    methods (Access=public)
        % see osCompute for details
        function obj = compute(obj, varargin)
            obj = osCompute(obj, varargin{:});
        end

        % see osPlot for details
        function uData = plot(obj, sensor, varargin)
            uData = osPlot(obj, sensor, varargin{:});
        end

        function val = timeAxis(obj)
            % The temporal samples for the lms filters
            val = ((1:size(linearFilters(obj.os, obj), 1)) - 1) ...
                * obj.timeStep;
        end

    end
end
