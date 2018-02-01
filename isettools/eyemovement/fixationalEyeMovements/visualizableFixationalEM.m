classdef visualizableFixationalEM < fixationalEM
    methods
        % Constructor
        function obj = visualizableFixationalEM(varargin)
            obj = obj@fixationalEM(varargin);
        end
        
        % Getter of propected properties (used by the demo app)
        function v = getValueOfProperty(obj, propertyName)
            v = obj.(propertyName);
        end
    end
end

