classdef visualizableFixationalEM < fixationalEM
% Create a visualizable fixational eye movement object
%
% Syntax:
%   visFixEM = visualizableFixationalEM();
%
% Description:
%    Create an instance of a visualizable fixational eye movement object.
%
% Inputs:
%    None required.
%
% Outputs:
%    The created visualizable fixational eye movement object.
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% See Also:
%    fixationalEM
%

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
            %    See pairs in header above.
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
            %    v            - VARIES. The property value.
            %
            % Optional key/value pairs:
            %    None.
            %

            v = obj.(propertyName);
        end
    end
end
