classdef osDisplayRGB < outerSegment 
% displayRGB subclass of the outersegment object
%
% Syntax:
%   os = osDisplayRGB();
%
% Description:
%    Bypasses the temporal filtering of the other outer segment subclasses
%    and passes the cone isomerizations without modification. It is
%    intended to be used for stimulus-referred retinal ganglion cell models
%    initially.
%
%    Examples are contained in the code. To access, type 'edit
%    osDisplayRGB.m' into the Command Window.
% 
% Inputs:
%    None required.
%
% Outputs:
%    The created outersegment object
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    subclasses: osLinear.m, osBioPhys.m
%

% History:
%    xx/xx/16  JRG/HJ/BW  ISETBIO Team, 2016
%    02/14/18  jnm        Formatting

% Examples:
%{
    displayRGBOS = osDisplayRGB();
%}

    % Public properties.
    properties (SetAccess = public, GetAccess = public)
    end

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        %rgbData - stored RGB data
        rgbData
    end

    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end

    % Public methods
    methods

        % Constructor
        function obj = osDisplayRGB(varargin)
            % Initialize the parent class
            %
            % Syntax:
            %   obj = osDisplayRGB([varargin])
            %
            % Description:
            %    Initialize the parent function.
            %
            % Inputs:
            %    varargin - (Optional) Any additional parameters provided.
            %
            % Outputs:
            %    obj      - The created outersegment Display RGB object.
            %
            % Optional key/value pairs:
            %    None.
            %
            obj = obj@outerSegment();

            % Initialize ourselves
            obj.rgbData = [];
        end

        % set function, see osIdentitySet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end

        % get function, see osIdentityGet for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end

    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sceneRGB)
            obj = osCompute(obj, sceneRGB); 
        end

        function plot(obj, sensor)
            osPlot(obj, sensor);
        end
    end

    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end

end
