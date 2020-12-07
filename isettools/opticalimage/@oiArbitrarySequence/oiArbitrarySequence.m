classdef oiArbitrarySequence < handle
% Class for generating a temporal sequence of optical images
%
% Syntax:
%   oiSequence = oiArbitrarySequence(oiList, oiTimeAxis, [varargin])
%
% Description:
%    Convenience class to be used similary to oiSequence for an arbitrary
%    collection of optical images
%
%    To create an oiArbitrarySequence using this function:
%
%        oiSeq = oiArbitrarySequence(oiList, timeAxis)
%
%    There are examples contained in the code. To access, type 'edit
%    oiSequence.m' into the Command Window.
%
% Inputs:
%    oiList             - Cell array of optical images
%    oiTimeAxis         - Vector. The time axis.
%
% Outputs:
%    The created OI Sequence object.
%
% Optional key/value pairs:
%    None
%
% See Also:
%   oisCreate, t_oiSequence, t_simplePhotocurrentComputation
%

% History:
%    06/28/2020  NPC  ISETBIO TEAM, 2020

properties

end

properties (Dependent)
    % length - Length of the sequence
    length
end

properties (SetAccess = private)
    % oiList - the cell array of optical images
    oiList
    
    % timeAxis - the oiSequence timebase
    timeAxis

end

methods  % public methods
    % constructor
    function obj = oiArbitrarySequence(oiList, oiTimeAxis, varargin)
        % Initialize the arbitrary optical image  object
        %
        % Syntax:
        %   obj = oiSequence(oiFixed, oiModulated, oiTimeAxis,
        %       modulationFunction, [varargin])
        %
        % Description:
        %    The constructor for the oi arbitrary sequence object.
        %
        % Inputs:
        %    See the class header.
        %
        % Outputs:
        %    obj                - Object. The created oiArbitrarySequence.
        %
        % Optional key/value pairs:
        %    None
        
        p = inputParser;
        p.addRequired('oiList', @iscell);
        p.addRequired('oiTimeAxis', @isnumeric);
        p.parse(oiList, oiTimeAxis, varargin{:});

        obj.oiList = p.Results.oiList;
        obj.timeAxis = p.Results.oiTimeAxis;

        assert(numel(obj.oiList) == numel(obj.timeAxis), 'time axis and oiList must have the same number of elements');
    end

    %% Define methods in the @oiSequence directory
    % Method to compute the maximum number of eye movement for current
    % sequence and a given integrationTime
    maxEyeMovementsNum = maxEyeMovementsNumGivenIntegrationTime(obj, ...
        integrationTime, varargin);

    % Method to return the oi at at the desired frame index
    oiFrame = frameAtIndex(obj, index);

    % Visualize the sequence
    [uData, hFig] = visualize(obj, varargin);

    % The time step of the timeAxis.
    function val = timeStep(obj)
        val = obj.timeAxis(2) - obj.timeAxis(1);
    end

    %% Local get methods
    % Return the length of the oiSequence
    function val = get.length(obj)
        val = numel(obj.timeAxis);
    end


    % Return the timeAxis used
    function val = get.timeAxis(obj)
        val = obj.timeAxis;
    end

end

end
