classdef rgcMosaic < cellMosaic
    %RGCMOSAIC - generates an rgcMosaic
    %
    %    rgcM = rgcMosaic(
    % The rgcMosaic class defines a particular RGC type and computational
    % model. The model specifies how to compute the RGC response from the
    % input. The type specifies cell and model parameters.
    %
    % The rgcLayer class  holds a collection of rgcMosaics.
    % The mosaics themselves are typically created by a function
    % of the rgcLayer class, such as
    %
    %   rgcLayer.mosaicCreate('model','LNP','type','cell type');
    %
    % Inputs:
    %    model: 'LNP','GLM' [subclasses of rgcMosaic]
    %    type: 'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget', 'Small Bistratified'
    %
    % Notes:
    %
    %   The RGC models are detailed in Chichilnisky & Kalmar, J. Neurosci (2002);
    %   Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
    %   and Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
    %   Nature (2008).
    %
    %   The computational model implemented for the coupled GLM model relies on
    %   code by <http://pillowlab.princeton.edu/code_GLM.html Pillow>, which is
    %   distributed under the GNU General Public License.
    %
    % See also: rgcLNP.m, rgcGLM.m, irCreate, s_initRetina
    %
    % Example:
    %
    %   ir.mosaicCreate('model','LNP','type','on midget');
    %
    %  ISETBIO wiki: <a href="matlab:
    %  web('https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cells','-browser')">RGCS</a>.
    %
    % JRG/BW ISETBIO team, 2015
    
    
    %% Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
    
    % Was Protected properties.  Changing to Public for debugging, and
    % maybe forever
    properties (SetAccess = public, GetAccess = public)
        
        % RFDIAMETER receptive field center diameter in MICRONS
        rfDiameter = [];
        
        %TCENTER temporal impulse response of the center in dt steps or 1
        %ms??
        tCenter =[];
        
        %TSURROUND  and of the surround (1 ms timing by default)
        tSurround = [];
        
        %TONICDRIVE baseline term for linear response; if nonzero, cell
        %spikes with no input
        tonicDrive;
        
        %RFDIAMAGNITUDE for making movies of response
        % rfDiaMagnitude;
        
        %RESPONSELINEAR Store the linear response after convolution
        responseLinear = [];
        
        %RESPONSESPIKES Store the spike times of the responses
        responseSpikes = [];
        
        %ELLIPSEMATRIX Store the parameters for the RGC sRF ellipses
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
        function obj = rgcMosaic(rgcLayer, cellType, varargin)
            %% Initialize an rgcMosaic for a particular cell type
            %
            %  rgcMosaic(rgcLayer, cellType, 'input mosaic',val)
            %
            % rgcMosaic is intialized based on a cell type and the
            % properties of the bipolar mosaic input. The bipolar mosaic is
            % specified by its index in the bpLayer object that is attached
            % to rgcLayer.input.
            %
            % BW, ISETBIO Team, 2017
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addRequired('rgcLayer',@(x)(isequal(class(x),'rgcLayer')));
            p.addRequired('cellType',@(x)(ismember(ieParamFormat(x),obj.validCellTypes)));
            p.addParameter('inMosaic',1,@isscalar);
            
            p.parse(rgcLayer,cellType,varargin{:});
            inMosaic = p.Results.inMosaic;
            
            % We need the parameters in the layer often enough.
            obj.parent = rgcLayer;
            
            % This the rgc mosaic type
            obj.cellType = strrep(lower(cellType),' ','');
            
            % Generate spatial RFs of the appropriate size for the cell type and TEE
            % Sets sRFcenter, sRFsurround
            obj.initSpace(cellType, 'inMosaic', inMosaic, varargin{:}); 
            
            % Sets temporal RF properties of tCenter/tSurround
            obj.initTime(rgcLayer);
            

        end
        
        % set function, see mosaicSet for details
        function obj = set(obj, varargin)
            mosaicSet(obj, varargin{:});
        end
        
        % get function, see mosaicGet for details
        function val = get(obj, varargin)
            val = mosaicGet(obj, varargin{:});
        end
                
    end
    
    properties (Constant)
        % VALIDCELLTYPES
        validCellTypes = {'onparasol','offparasol','onmidget','offmidget','onsbc'};
    end
    
    %% Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function window(obj)
            % Tip: Retrieve guidata using
            %    gui = guidata(obj.figureHandle);
            %
            obj.fig = rgcMosaicWindow(obj);
            
        end
        
        function val = timeAxis(obj)
            % Time steps in seconds.  Usually, dt is in 0.1 ms
            val = obj.dt*(1:length(obj.tCenter{1}))*1e-3;
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
            str = sprintf('Model: %s\n',class(obj));
            txt = sprintf('Cell type: %s\n',obj.cellType);
            str = addText(str,txt);
            
            % Mosaic properties
            txt = sprintf('N Trials %d\n',parent.nTrials);
            str = addText(str,txt);
            txt = sprintf('Patch size %d (um)\n',1e6*parent.size);
            str = addText(str,txt);
            
            % Spatial temporal properties
            txt = sprintf('Row,Col: %d, %d\n', size(obj.cellLocation));
            str = addText(str,txt);
            txt = sprintf('Time samples: %d\n',size(obj.responseLinear,3));
            str = addText(str,txt);
            
            txt = sprintf('Duration: %.0f ms\n',1e3*obj.dt*size(obj.responseLinear,3));
            str = addText(str,txt);
            
        end
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
