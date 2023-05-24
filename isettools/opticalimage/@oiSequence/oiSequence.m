classdef oiSequence < handle
% Class for generating a temporal sequence of optical images
%
% Syntax:
%   oiSequence = oiSequence(oiFixed, oiModulated, oiTimeAxis, ...
%       modulationFunction, [varargin])
%
% Description:
%    The stimuli for many psychophysical and physiological experiments can
%    be described as a fixed background image combined with a time-varying
%    image. The oiSequence class includes methods that combine the such
%    image combinations across time.
%
%    There is a companion method, oisCreate() that produces certain classic
%    oiSequences (harmonics, vernier targets, impulses). For some cases,
%    that function may provide a simpler interface for creating an
%    oiSequence.
%
%    To create an oiSequence using this function:
%
%        oiSeq = oiSequence(oiBackground, oiModulated, timeAxis, ...
%                           weightsTime, [varargin])
%
%    There are examples contained in the code. To access, type 'edit
%    oiSequence.m' into the Command Window.
%
% Inputs:
%    oiFixed            - Struct. The fixed optical image (background)
%    oiModulated        - Struct. A uniform optical image field.
%    oiTimeAxis         - Vector. The time axis.
%    modulationFunction - Vector. The modulation function.
%
% Outputs:
%    The created OI Sequence object.
%
% Optional key/value pairs:
%    composition        - String. Select one of the appropriate strings to
%                         indicate how to handle the fixed & modulated OIs.
%                         The options are: 'add', 'blend', and 'xor'.
%    modulationRegion   - Struct. A select region of interest to be used
%                         for the modulation.
%
% Notes:
%    * TODO: (Programming) Maybe set up some oiGet/Set routines with
%      the syntax such as the following:
%        oiSequence.get('oiFixed mumble') and
%        oiSequence.get('oiModulated mumble')?
%        Maybe oiSequence.get('frame', val)?
%
% See Also:
%   oisCreate, t_oiSequence, t_simplePhotocurrentComputation
%

% History:
%    xx/xx/16  NPC  ISETBIO TEAM, 2016
%    03/27/18  jnm  Formatting
%    07/02/19  JNM  Formatting update

% Examples:
%{
    % Harmonic
    oi = oiCreate('wvf human');
    params.freq = 6;
    params.contrast = 0.6;
    scene = sceneCreate('harmonic', params);
    oiModulated =  oiCompute(oi, scene);

    params.contrast = 0;  % contrast of the two frequencies
    scene = sceneCreate('harmonic', params);
    oiBackground =  oiCompute(oi, scene);
    stimWeights = fspecial('gaussian', [1, 50], 15);
    sampleTimes = 0.002 * (1:length(stimWeights));
    oiHarmonicSeq = oiSequence(oiBackground, oiModulated, ...
        sampleTimes, stimWeights/max(stimWeights), 'composition', 'blend');
    oiHarmonicSeq.visualize('movie illuminance');
 %}

properties
    % photonsFixed - Photons stored as a double
    photonsFixed;

    % photonsModulated - Photons stored as a double.
    photonsModulated;
end

properties (Dependent)
    % length - Length of the sequence
    length
end

properties (SetAccess = private)
    % oiFixed - the fixed oi (an oi, the background)
    oiFixed

    % oiModulated - the modulated oi (an oi, the pedestal)
    oiModulated

    % timeAxis - the oiSequence timebase
    timeAxis

    % Composition - Add or blend oiModulated to oiFixed?
    %   Whether to add the oiModulated to the oiFixed or to blend it with
    %   the oiFixed.
    composition;

    % modulationFunction - The modulating function
    %   This is an array of modulation values, one for each frame.
    modulationFunction;

    % modulationregion - The modulating region
    %   This is a struct describing the region extent to be modulated
    %   (for now just a radius).
    modulationRegion;
end

methods  % public methods
    % constructor
    function obj = oiSequence(oiFixed, oiModulated, oiTimeAxis, ...
            modulationFunction, varargin)
        % Initialize the optical image temporal sequence object
        %
        % Syntax:
        %   obj = oiSequence(oiFixed, oiModulated, oiTimeAxis,
        %       modulationFunction, [varargin])
        %
        % Description:
        %    The constructor for the oi sequence object.
        %
        % Inputs:
        %    See the class header.
        %
        % Outputs:
        %    obj                - Object. The created oiSequence.
        %
        % Optional key/value pairs:
        %    modulationRegion   - Struct. A structure describing the region
        %                         to be modulated.
        %    composition        - String. One of 'add', 'blend', or 'xor'.
        %                         Default 'add'.
        %    **Remainder needs to be filled out**
        %
        defaultModulationRegion = struct(...
            'radiusInMicrons', nan);

        p = inputParser;
        p.addRequired('oiFixed', @isstruct);
        p.addRequired('oiModulated', @isstruct);
        p.addRequired('oiTimeAxis', @isnumeric);
        p.addRequired('modulationFunction', @isnumeric);
        p.addParameter('modulationRegion', ...
            defaultModulationRegion, @isstruct);
        p.addParameter('composition', 'add', ...
            @(x)ismember(x, {'add', 'blend', 'xor'}));
        p.parse(oiFixed, oiModulated, oiTimeAxis, ...
            modulationFunction, varargin{:});

        obj.oiFixed = p.Results.oiFixed;
        obj.oiModulated = p.Results.oiModulated;
        obj.timeAxis = p.Results.oiTimeAxis;
        obj.modulationFunction = p.Results.modulationFunction;
        obj.modulationRegion = p.Results.modulationRegion;
        obj.composition = p.Results.composition;

        % The timeAxis must be the same length as the number of values in
        % the modulationFunction. If it only has 1 value, we are going to
        % assume that the value is delta T and we will create the whole
        % vector. If it has a vector, but that vector is not the same
        % length as the modulationFunction, we throw an error.
        if length(obj.timeAxis) == 1
            % First moment in time is 0. Increments by the set value.
            obj.timeAxis = obj.timeAxis * ...
                (0:(length(modulationFunction)) - 1);
        elseif length(obj.timeAxis) ~= length(obj.modulationFunction)
            error('Time axis does not match modulation function');
        end

        % Make sure that oiFixed and oiModulated have identical shape
        oiFixedSpatialSupport = round(oiGet(obj.oiFixed, ...
            'spatial support', 'microns'), 7);
        oiModulatedSpatialSupport = round(oiGet(obj.oiModulated, ...
            'spatial support', 'microns'), 7);

        if (any(size(oiFixedSpatialSupport) ~= ...
                size(oiModulatedSpatialSupport)))
            error(['Mismatch between spatial dimensions of ' ...
                'oiFixed, oiModulated']);
        end
        if (any(oiFixedSpatialSupport(:) ~= ...
                oiModulatedSpatialSupport(:)))
            error(['Mismatch between spatial support of ' ...
                'oiFixed, oiModulated']);
        end
    end

    %% Define methods in the @oiSequence directory
    % Method to compute the maximum number of eye movement for current
    % sequence and a given integrationTime
    maxEyeMovementsNum = maxEyeMovementsNumGivenIntegrationTime(obj, ...
        integrationTime, varargin);

    % Method for on-the-fly computation of the oi at desired index
    oiFrame = frameAtIndex(obj, index);

    % Visualize the sequence
    [uData, hFig] = visualize(obj, varargin);

    % The breadth of the timeAxis.
    function val = timeStep(obj)
        val = obj.timeAxis(2) - obj.timeAxis(1);
    end

    %% Local get methods
    % Return the length of the oiSequence
    function val = get.length(obj)
        val = numel(obj.modulationFunction);
    end

    % Return the modulationFunction used
    function val = get.modulationFunction(obj)
        val = obj.modulationFunction;
    end

    % Return the timeAxis used
    function val = get.timeAxis(obj)
        val = obj.timeAxis;
    end

    % Return the composition type used
    function val = get.composition(obj)
        val = obj.composition;
    end

end

end
