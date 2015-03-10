function optics = opticsCreate(opticsType,varargin)
% Create an optics structure
%
%   optics = opticsCreate(opticsType,varargin)
%
% The optics structure contains a variety of parameters,
% such as f-number  and focal length.
%
% Optics structures do not contain a spectrum structure.  Rather this is
% stored in the optical image that also holds the optics information.
%
% For diffraction-limited optics, the only parameter that matters really is
% the f-number.  The names of the standard types end up producing a variety
% of sizes that are only loosely connected to the names.
%
%      {'default', 'standard (1/4-inch)'}
%      {'standard (1/3-inch)'}
%      {'standard (1/2-inch)'}
%      {'standard (2/3-inch)'}
%      {'standard (1-inch)'}
%         
% There is one special case, human optics.  This creates an optics
% structure with human OTF data.
%
%      {'human'}            % Uses Marimont and Wandell (Hopkins) method
%      {'Ijspeert'}         % Ijspeert OTF (Not yet implemented)
%
% Example:
%   optics = opticsCreate('standard (1/4-inch)');
%   optics = opticsCreate('standard (1-inch)');
%
%   optics = opticsCreate('human');        % 3mm diameter is default
%   optics = opticsCreate('human',0.002);  % 4 mm diameter
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('opticsType'), opticsType = 'default'; end

opticsType = ieParamFormat(opticsType);

switch lower(opticsType)
    case {'default', 'standard(1/4-inch)', 'quarterinch'}
        optics = opticsQuarterInch;
    case 'human'
        % Pupil radius in meters.  Default is 3 mm
        if ~isempty(varargin), pupilRadius = varargin{1}; 
        else                   pupilRadius = 0.0015; % 3mm diameter default
        end
        % This creates a shift-invariant optics. The other standard forms
        % are diffraction limited.
        optics = opticsHuman(pupilRadius);
        optics = opticsSet(optics, 'model', 'shiftInvariant');
        optics = opticsSet(optics, 'name', 'human-MW');

    case 'mouse'
        % Pupil radius in meters.  
        % Dilated pupil : 1.009mm = 0.001009m
        % Contracted pupil : 0.178 mm
        % (Source : From candelas to photoisomerizations in the mouse eye by 
        % rhodopsin bleaching in situ and the light-rearing dependence of 
        % the major components of the mouse ERG, Pugh, 2004)
        % We use a default value, in between : 0.59 mm.
        if ~isempty(varargin)
            pupilRadius = varargin{1};
            if pupilRadius > 0.001009 || pupilRadius < 0.000178
                warning('Poor pupil size for the  mouse eye.')
            end
        else
            pupilRadius = 0.00059;   % default : 0.59 mm
        end
        % This creates a shift-invariant optics.  The other standard forms
        % are diffraction limited.
        optics = opticsMouse(pupilRadius);
        optics = opticsSet(optics,'model','shiftInvariant');
        
    case 'ijspeert'
        disp('Ijspeert optics not yet implemented')
        return;
    case {'standard(1/3-inch)','thirdinch'}
        optics = opticsThirdInch;
    case {'standard(1/2-inch)','halfinch'}
        optics = opticsHalfInch;
    case {'standard(2/3-inch)','twothirdinch'}
        optics = opticsTwoThirdInch;
    case {'standard(1-inch)','oneinch'}
        optics = opticsOneInch;
    otherwise
        error('Unknown optics type.');
end

% Default computational settings for the optical image
optics = opticsSet(optics,'offAxisMethod','cos4th');
optics.vignetting =    0;   % Pixel vignetting is off

end

%---------------------------------------
function optics = opticsQuarterInch
% Standard optics have a 46-deg field of view degrees

optics.type = 'optics';
optics = opticsSet(optics,'name','standard (1/4-inch)');
optics = opticsSet(optics,'model','diffractionLimited');

% Standard 1/4-inch sensor parameters
sensorDiagonal = 0.004;
FOV = 46;
fLength = inv(tan(FOV/180*pi)/2/sensorDiagonal)/2;

optics = opticsSet(optics,'fnumber',4);  % focal length / diameter
optics = opticsSet(optics,'focalLength', fLength);  
optics = opticsSet(optics,'otfMethod','dlmtf');

end

%---------------------------------------
function optics = opticsThirdInch
% Standard 1/3-inch sensor has a diagonal of 6 mm
%
optics.type = 'optics';
optics = opticsSet(optics,'name','standard (1/3-inch)');
optics = opticsSet(optics,'model','diffractionLimited');

optics = opticsSet(optics,'fnumber',4);

% Standard optics have a 46-deg field of view degrees
FOV = 46;
sensorDiagonal = 0.006;
fLength = inv(tan(FOV/180*pi)/2/sensorDiagonal)/2;

optics = opticsSet(optics,'focalLength', fLength);  
optics = opticsSet(optics,'otfMethod','dlmtf');

end

%---------------------------------------
function optics = opticsHalfInch
optics.type = 'optics';
optics = opticsSet(optics,'name','standard (1/2-inch)');
optics = opticsSet(optics,'model','diffractionLimited');

optics = opticsSet(optics,'fnumber',4);  % focal length / diameter

% Standard optics have a 46-deg field of view degrees
FOV = 46;
sensorDiagonal = 0.008;
fLength = inv(tand(FOV)/2/sensorDiagonal)/2;

% Standard 1/2-inch sensor has a diagonal of 8 mm
optics = opticsSet(optics,'focalLength', fLength);  
optics = opticsSet(optics,'otfMethod','dlmtf');

end

%---------------------------------------
function optics = opticsTwoThirdInch
%

optics.type = 'optics';
optics = opticsSet(optics,'name','standard (2/3-inch)');
optics = opticsSet(optics,'model','diffractionLimited');

FOV = 46;
sensorDiagonal = 0.011;
fLength = inv(tan(FOV/180*pi)/2/sensorDiagonal)/2;

optics = opticsSet(optics,'fnumber',4);
optics = opticsSet(optics,'focalLength', fLength);  
optics = opticsSet(optics,'otfMethod','dlmtf');

end

%---------------------------------------
function optics = opticsOneInch
% Standard 1-inch sensor has a diagonal of 16 mm

optics.type = 'optics';
optics = opticsSet(optics,'name','standard (1-inch)');
optics = opticsSet(optics,'model','diffractionLimited');

FOV = 46;
sensorDiagonal = 0.016;
fLength = inv(tan(FOV/180*pi)/2/sensorDiagonal)/2;

optics = opticsSet(optics,'fnumber',4);
optics = opticsSet(optics,'focalLength', fLength);  
optics = opticsSet(optics,'otfMethod','dlmtf');
        
end

%---------------------------------------
function optics = opticsHuman(pupilRadius)
% We use the shift-invariant method for the human and add the OTF
% data to the OTF fields.   We return the units in cyc/mm.  We use 300
% microns/deg as the conversion factor.
% EC - 300um/deg corresponds to a distance of 17mm (human focal length)

% We place fnumber and focal length values that
% are approximate for diffraction-limited in those fields, too.  But they
% are not a good description, just the DL bounds for this type of a system.
%
% The pupilRadius should be specified in meters
%

if notDefined('pupilRadius'), pupilRadius = 0.0015; end
fLength = 0.017;  %Human focal length is 17 mm

optics.type = 'optics';
optics.name = 'human';
optics      = opticsSet(optics, 'model', 'shiftInvariant');

% Ratio of focal length to diameter.  
optics = opticsSet(optics,'fnumber',fLength/(2*pupilRadius));  
optics = opticsSet(optics,'focalLength', fLength);  

optics = opticsSet(optics, 'otfMethod', 'humanOTF');

% Compute the OTF and store it.  We use a default pupil radius, dioptric
% power, and so forth.

dioptricPower = 1/fLength;      % About 60 diopters

% We used to assign the same wave as in the current scene to optics, if the
% wave was not yet assigned.  
wave = opticsGet(optics, 'wave');

% The human optics are an SI case, and we store the OTF at this point.  
[OTF2D, frequencySupport] = humanOTF(pupilRadius, dioptricPower, [], wave);
optics = opticsSet(optics,'otfData',OTF2D);

% Support is returned in cyc/deg.  At the human retina, 1 deg is about 300
% microns, so there are about 3 cyc/mm.  To convert from cyc/deg to cyc/mm
% we divide by 0.3. That is:
%  (cyc/deg * (1/mm/deg)) cyc/mm.  1/mm/deg = 1/.3
frequencySupport = frequencySupport * (1/0.3);  % Convert to cyc/mm

fx     = frequencySupport(1,:,1);
fy     = frequencySupport(:,1,2);
optics = opticsSet(optics, 'otffx', fx(:)');
optics = opticsSet(optics, 'otffy', fy(:)');

optics = opticsSet(optics, 'otfWave', wave);

% figure; 
% mesh(frequencySupport(:,:,1),frequencySupport(:,:,2),OTF2D(:,:,20));
% mesh(abs(otf2psf(OTF2D(:,:,15))))
%
end