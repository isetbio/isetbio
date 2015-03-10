function res = rgcGet(rgcP,fieldName,varargin)
% Implements the get function for the independent rgcParameters class
%
%    res = rgcGet(rgcP,fieldName,[optParam])
%
% rgcP:     the rgc parameter object 
% fieldName: string describing the field that we want
%
% This routine should not be called directly. It is called from the
% class get function in the form: 
%
%      res = obj.get(fieldName,optParam);
%
% This routine also manages the interface to the get function for the
% layers structure.
%
% The list of rgcParameters and rgcLayers is below.  These may be moved
% from here to the rgcParameters file.
%
% RGC object parameters:
%     'name'  - Name of this particular instance of the rgcParameter class
%     'noise' - Noise function that uses meanV and stdV as parms
%     'meanv'    - noise function mean
%     'stdv'     - noise function standard deviation
%     'noiseframes' - Warm up noise frames?
%
%     'data' - Absorptions across the array and across time
%
%     'dt'   - Time step in milliseconds
%     'dT_default'  - If no sensor, this is the default
%     'feedbacktimes' - Temporal samples for feedback
%     'trdur' - Duration (ms) of stored impulse responses
%     'tshift' - Broken ...
%     'linf'   - Broken
%     'cptr'   - Coupling impulse response function - Broken
%     'conespacing' -        Stored in sensor units are microns
%     'conespacingdefault' - Only used when there is no sensor ...
%     'savedir'  - Broken?
%     'cutoff'   - Broken?
%     'distancefunction' - How to compute distance between RGCs
%     'layer'    -  Layer class
%     'nframe'   -  Number of time frames
%     'conegrid' -  Indices, without units, of cone grid.  Useless?
%     'nlayers'  -  Number of ganglion cell layers
%     'rf' - Receptive field?  Why here and not in layer instead?
%
%  Eye movement related - Removed (BW).
%     'eyemovetype' - Eye movement function
%     'maxeyemove'  - Upper bound
%     'eyemovestd'  - Std dev. of Brownian eye movement
%     'noeyemove'   - Turn off eye movements
%     'eyepath'     - Stored path?
%
% RGC Layer object parameters  - These are a mess as to what is returned.
%     'name'   - Name of this layer
%     'cellspacing' - Cell spacing in microns re: cone mosaic?
%     'gridsize'    - Cone grid size?
%     'cellgrid'    - p=rgcP.get{''}; foo = p{1}(:,1); foo{1}  (um ?)
%     'cellloc'     - p=rgcP.get('cellloc'); p{1}
%     'layercenter' - p=rgcP.get('layercenter'); p{1}
%     'parent'      - rgcP structure
%     'rfcov'
%     'rfstd'
%     'rfscoeffs'
%     'coneweights'
%     'trdur'
%     'centertr'
%     'centerbdiv'
%     'centerf'
%     'centernormf'
%     'centergamma'
%     'surroundtr'
%     'surroundbdiv'
%     'surroundf'
%     'surroundnormf'
%     'surroundgamma'
%     'gammaftr'
%     'gammabdiv'
%     'gammaf'
%     'gammanormf'
%     'gammagamma'
%     'dt'
%     'tshift'
%     'linf'
%     'rf'   - Receptive field
%     'actthreshparam' - Activity threshold parameter  (mv? or unscaled?)
%     'conespacing'    - In units of microns?  Meters?
%     'cutoff'
%     'wrange'         - Range of some voltage
%     'nrange'
%     'rgcvtimeseries'
%     'currentspikeintegral'
%     'currentspkts'
%     'currentlints'
%     'currentconnec'
%     'bordersize'
%     'overridesizedefault'
%     'overridensize'
%
% See also: rgcMapParameterField, rgcParameters, rgcLayers
%
% Examples:
%    rgcP = rgcParameters;
%    rgcP.addOnLayer();
%    rgcP.set('data',rand(50,50,8));
%
%  the preferred calls:
%
%    absorptions = rgcP.get('absorptions');
%    cVolts      = rgcP.get('cone voltages');
%    layer       = rgcP.get('layer',1);
%
% are equivalent to these calls
%
%    absorptions = rgcGet(rgcP,'cone absorptions');
%    cVolts      = rgcGet(rgcP,'cone voltages');
%    layer       = rgcGet(rgcP,'layer',1);
%
%
% If you want to add an alias to an rgc parameters field, edit the function
% rgcMapParameterField. Do not add aliases into the case statements in this
% function.
%
% (c) Stanford Vista Team, 2009

% TODO:  No error is returned when the parameter is not found.

if notDefined('rgcP'),  error('Parameter object needed'); end
if notDefined('fieldName'),  error('Field name needed'); end

% Map the input field name to a standard format, lower case and no spaces
fieldName = rgcMapParameterField(fieldName);

% correctfield is a standard way of calling the field
switch(fieldName)
    
    % Book-keeping parameters
    case {'name'}
        res = rgcP.name;
    case{'type'}
        % If this is the rgc object, class is 'rgcParameters'
        res = class(rgcP);
        
        % Not needed ... we should handle noise separately, I think.
    case {'noise'}
        res = rgcP.noise;
    case {'noiseframes'}
        res = rgcP.noise.nFrame;
        
        % ISET structures and cone input related
    case {'absorptions'}
        res = rgcP.absorptions;
        
    case {'conegrid'}
        % res = rgcP.get('coneGrid','mm')
        % Default is in microns
        if isempty(varargin), units = 'um';
        else units = varargin{1};
        end
        cSpacing = rgcP.get('cone spacing',units);
        gs       = rgcP.get('cone grid size');
        if isempty(gs), warning('No cone grid'), res = []; return; end

        cg{1} = (1:gs(1))*cSpacing; cg{1} = cg{1} - mean(cg{1});
        cg{2} = (1:gs(2))*cSpacing; cg{2} = cg{2}  -mean(cg{2});
        res = cg;

    case{'conegridsize'}
        % The row and column dimensions of the cones
        % rgcP.get('cone grid size')
        sensor = rgcP.sensor;
        if isempty(sensor), warning('No sensor'); res=[]; return; end
        res = sensorGet(sensor,'size');
        
    case {'conespacing'}
        % You can request the units.  Derived from the sensor structure.
        % cs = rgcP.get('cone spacing','um');
        % res = rgcP.coneSpacing;
        sensor = rgcP.get('sensor');
        pixel  = sensorGet(sensor,'pixel');
        if isempty(varargin), units = 'um';
        else units = varargin{1};
        end
        res = pixelGet(pixel,'height',units);
        
    case {'coneimagesize'}
        if isempty(varargin), units = 'um';
        else units = varargin{1};
        end
        res = rgcP.get('cone grid size')*rgcP.get('cone spacing',units);
        
    case {'cvolts'}
        % Cone volts?  Obsolete?
        res = rgcP.cVolts;
        
    case {'sensor'}
        res = rgcP.sensor;
        
    case {'oi'}
        res = rgcP.oi;
        
        
        
        % Timing parameters
    case {'dt'}
        % Returned by default in milliseconds
        % rgcP.get('dT')
        if isempty(varargin), units = 'ms'; 
        else units = varargin{1};
        end
        
        % Get timing from sensor, if available.
        sensor = rgcP.sensor;
        if isempty(sensor), warning('No sensor'); res = []; return; end
        res = sensorGet(sensor,'expTime',units);
               
    case{'trsamples'}
        % rgcP.get('trsamples') 
        % rgcP.get('trsamples','sec')
        % The set of temporal response samples (ms default) used for many
        % temporal coupling functions.

        if isempty(varargin), units = 'ms'; 
        else units = varargin{1};
        end
        dT  = rgcP.get('dT',units);
        res = (dT : dT : rgcP.get('trDur',units));

    case {'trdur'}
        % Temporal response function support (total duration)
        if isempty(varargin), units = 'ms';
        else units = varargin{1};
        end
        res = rgcP.trDur;
        % Stored in ms.  Correct to secs and then scale
        res = res * 10^-3 * ieUnitScaleFactor(units);
        
    case {'nframe'}
        % Number of temporal samples
        %         res = rgcP.nFrame;
        res = size(rgcP.get('cVolts'),3);

    case{'tshift'}
        % No idea.
        res = rgcP.tShift;
        
        % No idea
    case{'linf'}
        res = rgcP.linF;
         
        % Should go away, I think.  The noise process should be established
        % separately.
        %case {'meanv'}
        %         % The mean of what?
        %         res = rgcP.meanV;
        %
        %     % case {'stdv'}
        %         res = rgcP.stdV;

        % Probably should go away
    case {'savedir'}
        res = rgcP.saveDir;
        
        % Connections related?
    case {'cutoff'}
        res = rgcP.cutoff;    
    case {'distancefunction'}
        res = rgcP.distanceFunction;
    
        % Layers are a big part and separate class.
    case{'layer'}
        % rgcP.get('layer',2)
        % If layer is not specified, defaults to 1.
        if isempty(varargin)
            if (rgcP.get('nLayers') == 1), layerNbr = 1;
            else error('Specify layer number.')
            end
        else layerNbr = varargin{1};
        end
        res = rgcP.layers{layerNbr};
              
    case {'alllayers','layers'}
        % All are in a cell array, normally
        res = rgcP.layers;        
    case {'nlayers'}
        res = length(rgcP.layers);
        
        % Cone input grid related       
    case {'rf'}
        % This is should deleted.
        % res = rgcP.computeRF();
        error('Get an RF from the layer, not the rgcP');

    otherwise
        try
            % Maybe this is a layer attribute, building a cell
            nL = rgcP.get('nLayers');
            res = cell(1,nL);
            for ii = 1:nL
                res{ii} = rgcP.getLayer(ii).get(fieldName);
            end
            fprintf('field name %s is a layer parameter \n',fieldName);
        catch me 
            % No message in matlab 2008a
            warning(me.identifier,me.message)
        end

end
