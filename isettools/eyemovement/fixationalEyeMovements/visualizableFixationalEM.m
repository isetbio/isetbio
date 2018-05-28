classdef visualizableFixationalEM < fixationalEM
% Create a visualizable fixational eye movement object
%
% Syntax:
%   visFixEM = visualizableFixationalEM();
%
% Description:
%    Create an instance of a visualizable fixational eye movement object
%    which is a subclass of @fixationalEM. The only purpose of this class
%    is to gain access to protected properties of @fixationalEM which are
%    used to by the fixationalEMdemo.mlapp to visualize the internals of 
%    the em generation algorithm.
%
% Inputs:
%    None required.
%
% Outputs:
%    The created visualizable fixational eye movement object.
%
% Optional key/value pairs:
%    None at the moment.
%
% See Also:
%    fixationalEM
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments

    methods
        % Constructor
        function obj = visualizableFixationalEM(varargin)
            % Initialize the visualizableFixationalEM class
            %
            % Syntax:
            %   obj = visualizableFixationalEM([varargin])
            %
            % Description:
            %    The constructor for the visualizable fixational eye
            %    movement object.
            %
            % Inputs:
            %    None required.
            %
            % Outputs:
            %    obj - The visualizableFixationalEM object.
            %
            % Optional key/value pairs:
            %    None at the moment.
            %

            obj = obj@fixationalEM(varargin);
        end
     
        % Getter of propected properties (used by the demo app)
        function v = getValueOfProperty(obj, propertyName)
            % Property retriever, used by demo app
            %
            % Syntax:
            %   v = getValueOfProperty(obj, propertyName)
            %   v = obj.getValueOfProperty(propertyName)
            %
            % Description:
            %    This is a simple property retriever for the demo app.
            %
            % Inputs:
            %    Obj          - Object. The visualizableFixationalEM object
            %    propertyName - String. The property name
            %
            % Outputs:
            %    v            - The property value.
            %
            % Optional key/value pairs:
            %    None.
            %

            v = obj.(propertyName);
        end
    end
end
