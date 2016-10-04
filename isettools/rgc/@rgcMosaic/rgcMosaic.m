classdef rgcMosaic < handle
% Generates an rgcMosaic 
%
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
% JRG/BW ISETBIO team, 2015

    
    %% Define object
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        Parent;
    end
    
    % Protected properties.
    properties (SetAccess = protected, GetAccess = public)
        
        % The type of computational model for the RGC
        cellType;           % Possible types are listed in header
        rfDiameter;         % receptive field center diameter

        % We should estimate the rf center sigma
        % rfDiaMagnitude;
        
        % Cell array cellLocation{i}{j} = [x,y] position (microns)
        cellLocation;
        sRFcenter;          % spatial RF of the center on the receptor grid
        sRFsurround;        % spatial RF of the surround
        tCenter;            % temporal impulse response of the center
        tSurround;          %    and of the surround (1 ms timing by default)
        tonicDrive;         % DC term for linear response
        rfDiaMagnitude;     % for making movies of response
        responseLinear;     % Store the linear response after convolution
        responseSpikes;     % Store the spike times of the responses

    end
    
    properties (Access = public)
        % When we open the figure for the mosaic, we store the handle here
        figureHandle;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = rgcMosaic(ir, mosaicInd)
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
            obj.cellType = mosaicInd;
            switch ieParamFormat(mosaicInd)
                case{'onparasol'}
                    mosaicInd = 1;
                case{'offparasol'}
                    mosaicInd = 2;
                case{'onmidget'}
                    mosaicInd = 3;
                case{'offmidget'}
                    mosaicInd = 4;
                case{'smallbistratified','sbc'}
                    mosaicInd = 5;
                otherwise
                    error('Unknown cell type');
            end
            
            % Generate spatial RFs of the approrpiate size for the cell type and TEE
            obj.rgcInitSpace(ir, mosaicInd);
            obj.rgcInitTime(ir, mosaicInd);
            
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
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function window(obj)
            obj.figureHandle = mosaicWindow(obj);
            % Tip: Retrieve guidata using
            %    gui = guidata(obj.figureHandle);
            %
        end
        
        function str= describe(obj)
            % Print the relevant text to a string
            % This is used in the display window
            
            str =sprintf('Cell type: %s\n',obj.cellType);
            txt = sprintf('Model: %s\n',class(obj));
            str = addText(str,txt);
            txt = sprintf('N Trials %d\n',obj.Parent.numberTrials);
            str = addText(str,txt);
            
            parent = obj.Parent;
            txt = sprintf('Row,Col: %d, %d\n',parent.row,parent.col);
            str = addText(str,txt);
            txt = sprintf('Patch size %d (um)\n',1e6*parent.spacing);
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
