classdef rgcMosaic < handle
%RGCMOSAIC - generates an rgcMosaic 
% Each RGC mosaic has a particular model and a type.  The model specifies
% how we compute the RGC response, and the type specifies the parameters of
% the model given the type.
%
% The inner retina object holds a collection of retinal ganglion cell
% mosaics. The mosaics themselves are typically created by a function of
% the inner retina class, such as
%
%   ir.mosaicCreate('model','LNP','type','cell type'); 
%
% This file contains the constructor and the window call.  Other functions,
% such as plot(), are separated in the @rgcMosaic directory.
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
        Parent;
    end
    
    % Was Protected properties.  Changing to Public for debugging, and
    % maybe forever
    properties (SetAccess = public, GetAccess = public)
        
        %CELLTYPE The type of computational model for the RGC
        cellType = 'onparasol';           % Possible types are listed in header
        
        % RFDIAMETER receptive field center diameter in MICRONS
        rfDiameter = [];
               
        %CELLLOCATION Cell array cellLocation{i}{j} = [x,y] position (microns)
        cellLocation;
        
        %SRFCENTER spatial RF of the center on the cone mosaic grid
        sRFcenter = [];           
        
        %SRFSURROUND spatial RF of the surround
        sRFsurround = [];        
        
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
        %FIGUREHANDLE When we open the figure for the mosaic, we store the handle here
        figureHandle;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
        
    end
    
    %% Public methods
    methods
        
        % Constructor
        function obj = rgcMosaic(ir, mosaicInd, varargin)
            %% Initialize an rgcMosaic for a particular cell type
            %
            %       initialize(obj, innerRetina, cellType)
            %           [only called internally from rgcMosaic.m]
            %
            % The object is intialized based on a series of input parameters that
            % can include the location of the retinal patch.
            %
            % First, the name of the cell type is assigned based on the value passed in
            % the type parameter. Next, spatial receptive fields of the appropriate
            % size are generated for the array of RGCs of that particular type. Then
            % the RGB temporal impulse responses for the center and surround are
            % generated.
            %
            % Switch cell type string to index number
            % The index number helps with the generation of the receptive fields and
            % impulse responses of the appropriate parameters for the cell type.
            obj.cellType = strrep(lower(mosaicInd),' ','');
            
            % Generate spatial RFs of the appropriate size for the cell type and TEE
            obj.rgcInitSpace(ir, mosaicInd,varargin{:}); % Sets sRFcenter, sRFsurround
            
            % Sets temporal RF properties of tCenter/tSurround
            obj.rgcInitTime(ir);             
            
            % We need the parameters in the parent often enough.  So put in
            % a pointer to it here.
            obj.Parent = ir;
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
    
    %% Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function window(obj)
            obj.figureHandle = mosaicWindow(obj);
            % Tip: Retrieve guidata using
            %    gui = guidata(obj.figureHandle);
            %
        end
        
        function val = timeAxis(obj)
            % Time steps in seconds.  Usually, dt is in 0.1 ms
            val = obj.dt*(1:length(obj.tCenter{1}))*1e-3;
        end
        
        
        
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
