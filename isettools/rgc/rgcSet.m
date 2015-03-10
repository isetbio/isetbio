function rgcP = rgcSet(rgcP,fieldName,val,varargin)
%Implements the set function for the rgcParameters class
%
%   res = rgcSet(rgcP,fieldName,val)
%
%  rgcP:     the parameter object (from the class rgcParameters)
%  fieldName: string describing the field that we want
%  val:     val to which we are setting the parameter
%
% This function is equivalent to
%
%   obj = obj.set(fieldName,val);
%
% which is probably the way it should be called
% for the layers you have special functions
%
% Example:
%    rgcP = rgcParameters;
%    rgcSet(rgcP,'absorptions',rand(100,100,5));
%
%  is equivalent to
%    rgcP.set('absorptions',rand(100,100,5));
%
%  but to add layers, you have to use
%    rgcP.addLayer();
%
% If you feel like an alias is missing for a field, feel free to add it in
% rgcMapParameterField But not in this function, which only uses the
% standard name.
%
% (c) Stanford VISTA Team, 2010

if notDefined('rgcP'), error('rgcParameter object needed'); end
if notDefined('fieldName'), error('Field name needed');   end

% Map multiple possible field names to a single standard field
fieldName = rgcMapParameterField(stringFormat(fieldName));

% correctfield is the standard way of calling the field
switch(fieldName)
    case {'name'}
        rgcP.name = val;
        
        % Cone related
    case {'absorptions'}
        % Cone absorptions computed by ISET
        rgcP.absorptions = val;
        
    case {'cvolts'}
        % Volts, probably after adaptation?
        rgcP.cVolts = single(val);
        
    case {'conespacing'}
        % Separation between cones.  But I think this should really be in
        % sensor, not here. It should be a sensorGet when you ask for it,
        % and it should never be set here.
        rgcP.coneSpacing = val;
        
    case {'conespacingdefault'}
        % This should go away, too, as per above.
        rgcP.coneSpacingDefault = val;
        
        % ISET related
    case {'scene'}
        % ISET scene
        rgcP.scene = val;
        
    case {'sensor'}
        % ISET sensor
        rgcP.sensor = val;
        
    case {'oi'}
        % ISET optical image
        rgcP.oi = val;
        
        % Temporal parameters
    case {'dt'}
        % Time step.  I think in ms.
        rgcP.dT = val;
        
    case {'trdur'}
        % ? Total duration of temporal samples?
        rgcP.trDur = val;
    case {'nframe'}
        % Probably should be deleted; get should derive from cVolts
        rgcP.nFrame = val;
        
        % Connection related
    case {'distancefunction'}
        rgcP.distanceFunction = val;
        
        % Spatial properties of RGC grid
    case {'gridsize'}
        rgcP.gridSize = val;
        
    case {'imgsize'}
        % Not sure if this is cVolt size or rgcSize.  Needs to be
        % clarified.
        rgcP.imSize = val;
        
        % Layers
    case{'layer'}
        error('see rgcParameters to handle layers.')
        
    case {'nlayers'}
        rgcP.nLayers = val;
        
    otherwise
        error('Unknown rgcParameter field: %s.',fieldName);
        
end

return