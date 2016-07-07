classdef rgcLNP < rgcMosaic
% rgcMosaic cell type with a LNP (linear-nonlinear-Poisson) computational model
%
% This model is found in Chichilnisky & Kalmar, J. Neurosci (2002) and
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli, Nature
% (2008).
% 
%  The computational model implemented here relies on code by
% <http://pillowlab.princeton.edu/code_GLM.html Pillow>, which is
% distributed under the GNU General Public License.
%
% rgcLNP is a subclass of rgcMosaic. It is called when creating a new LNP
% model rgcMosaic for an inner retina object.  Typically we get here from
% the inner retina object via a call
%
%   ir.mosaicCreate('model','LNP','type','your type goes here')
% 
% See also:
%
% Example:
%   os = osCreate('identity');        % A pass through from the stimulus
%   ir = irCreate(os,'name','myRGC'); % An inner retina container
%   ir.mosaicCreate('model','LNP','type','on midget'); % This  mosaic
%
% 9/2015 JRG    

    % Public, read-only properties.
    properties (SetAccess = public, GetAccess = public)
        
        % Should this be in rgcMosaic or in each of the special types of
        % rgcMosaics? (BW)
        
        % Number of repeats
        numberTrials = 10;
    end
           
    % Protected properties.
    properties (SetAccess = private, GetAccess = public)
        % The linear voltage is computed from the properties of the parent
        % class, which are attached below (rgcMosaic)
        % This includes the first voltage response, which is called
        %
        % responseLinear
               
        % Pillow promotes the linear input voltage using a nonlinear
        % function that he calls the generator function.  By default this
        % is an exponential.
        generatorFunction;
        
        % Parameter to specify the time bins Pillow uses for coupling and
        % post spike filters (10 ms default)
        dt = 0.01;
        
        % The nonlinear voltage response after application of the generator
        % function and the spike coupling responses is represented here
        responseVoltage;
        
        % The spikes are computed from the response voltage
        responseSpikes;

        % These hold the parameters used in the computation.
        % This is the response after a spike
        postSpikeFilter;
              
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcLNP(rgc, mosaicType)
            % Constructs the rgcLNP object
            %
            %       mosaic = rgcLNP(ir,mosaicType);
            %
            % The mosaic object is initialized with the generatorFunction,
            % postSpikeFilter and couplingFilter properties.
            %
            % The parent class is rgcMosaic, which specify the linear
            % convolutional properties.
            %
            % Example:
            %
            % See also: rgcMosaicCreate
            %
            % (c) isetbio
            % 09/2015 JRG
            
            
            % Initialize the parent class            
            obj = obj@rgcMosaic(rgc, mosaicType);

            % Effect of a spike on output voltages
            obj.generatorFunction = @exp;
            
            % Post spike filter
            obj.postSpikeFilter = buildPostSpikeFilter(obj.dt);
            
        end
        
        % set function, see for details
        function obj = set(obj, varargin)
            % obj = set@rgcMosaic(obj, varargin);
            mosaicSet(obj, varargin{:});
        end
        
        % get function, see for details
        function val = get(obj, varargin)
           val = mosaicGet(obj, varargin{:});
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        
                % Overloaded display function to order properties as we want to
        function disp(obj)
            
            % Get the string from the builtin disp and then rearrange
            s1 =evalc('builtin(''disp'', obj'')');
            
            % Find newline character locations
            newlineChars = strfind(s1,sprintf('\n'));
            
            % Pull out lines with properties and values
            for nlInd = 2:length(newlineChars)
                propertystr{nlInd} = s1(newlineChars(nlInd-1)+1:newlineChars(nlInd)-1);
            end
            
            % Display the reordered property list
            % Display the first line with class type
            disp(sprintf('\n'));
            disp(s1(1:newlineChars(1)))
            % Order determined through trial and error
            propertyOrder = [9:18 3:4 7:8 5:6];
            for nlInd = propertyOrder
                disp(propertystr{nlInd});
            end
        end
        
        initialize(obj);
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
