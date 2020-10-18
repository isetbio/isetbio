function [oType, p] = ieParameterOtype(param)
% Determine object type from parameter name
%
% Syntax:
%   [oType, p] = ieParameterOtype(param, varargin)
%
% Description:
%    Parse the parameter string to determine which object (oType) and the
%    parameter itself (p).
%
% Inputs:
%    param - The input string which to parse
%
% Outputs:
%    oType - The object
%    p     - The parameter
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note:  BW - Removed  ISET specific cases.  ]
%    * [Note: JNM - rtpsfsize was listed 2x in the switch case for
%      ieParamFormat. I left it in oi and removed it from optics (based on
%      order of case appearance. If this is not correct, please let me know
%      and I will change which it has been left in.]
%    * TODO: Complete TODO listed below.
%

% History:
%    xx/xx/13       Copyright Imageval LLC, 2013
%    12/12/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    [oType, p] = ieParameterOtype('optics/fnumber')
%}
%{
    % This one won't work because size is available for too many objects
	ieParameterOtype('size')
%}
%{
  [o p] = ieParameterOtype('oi size')
  [o p] = ieParameterOtype('scene/hfov')
  [o p] = ieParameterOtype('scene hfov')
  [o p] = ieParameterOtype('optics/fnumber')
%}


%%
if ~exist('param', 'var') || isempty(param), error('Param required.'); end

%% If the param is one of the objects, return the object string as the oType

p = [];
switch ieParamFormat(param)
    case 'scene'
        oType = 'scene'; return;
    case 'oi'
        oType = 'oi'; return;
    case 'optics'
        oType = 'optics'; return;
    case 'lens'
        oType = 'lens'; return;
    case {'em', 'eyemove', 'eyemovement'}
        oType = 'em'; return;
end

%% Find the string before the first space or the first '/'. 

% Special characters are spaces and /
c1 = strfind(param, ' ');   % Find the spaces
c2 = strfind(param, '/');   % Find the '/'

% Find the first space or '/'. If we checked for whole string, we wouldn't
% need the bit above. 
pos = min([c1, c2]);

% Parse and return the string as oType
% TODO:  We will need to add 'cone' and 'macular'
oType = [];
if ~isempty(pos)
    switch param(1:(pos-1))
        case 'scene'
            oType = 'scene';
        case 'oi'
            oType = 'oi';
        case 'optics'
            oType = 'optics';
        case 'lens'
            oType = 'lens';
        % case 'sensor'
            % oType = 'sensor';
        % case 'pixel'
            % oType = 'pixel';
        % % case {'vci', 'ip'}
            % % oType = 'ip';
        case {'em', 'eyemove', 'eyemovement'}
            oType = 'em';
    end
    
    % Check for success. Return the parameter, without the prepended term
    if ~isempty(oType), p = param((pos + 1):end);  return; end
end

%% We didn't find the oType yet.

% Maybe param is one of the unique object parameter strings. We have a
% list here.

% By the way ...
% I think there is a better way to do this in vistasoft, as per Adrian's
% coding using hashing. Not sure we can use it in previous versions of
% Matlab, though.
%
% TODO:
% Need to decide where the rt/raytrace goes. Not organized correctly yet.
% A lot of these should go away because they are related to objects not in
% ISETBIO, but only in ISET.
p = param;
switch ieParamFormat(param)
    % TODO
    %     case {''}
    %         oType = 'lens';
    %     case {''}
    %         oType = 'macular';
    %     case {''}
    %         oType = 'cone';
    case {'objectdistance', 'meanluminance', 'luminance', ...
            'illuminant', 'illuminantname', 'illuminantenergy', ...
            'illuminantphotons', 'illuminantxyz', 'illuminantwave', ...
            'illuminantcomment', 'illuminantformat'}
        oType = 'scene';
        
    case {'optics', 'opticsmodel', 'diffusermethod', 'diffuserblur'...
            'psfstruct', 'sampledrtpsf', 'psfsampleangles', ...
            'psfanglestep', 'psfimageheights', 'raytraceopticsname', ...
            'rtpsfsize'}
        oType = 'oi';
        
    case {'fnumber', 'effectivefnumber', 'focallength', 'power', ...
            'imagedistance', 'imageheight', 'imagewidth', ...
            'numericalaperture', 'aperturedameter', 'apertureradius', ...
            'magnification', 'pupilmagnification', 'offaxismethod', ...
            'cos4thmethod', 'cos4thdata', 'otfdata', 'otfsize', ...
            'otffx', 'otffy', 'otfsupport', 'psfdata', 'psfspacing', ...
            'psfsupport', 'incoherentcutoffspatialfrequency', ...
            'maxincoherentcutoffspatialfrequency', 'rtname', ...
            'raytrace', 'rtopticsprogram', 'rtlensfile', ...
            'rteffectivefnumber', 'rtfnumber', 'rtmagnification', ...
            'rtreferencewavelength', 'rtobjectdistance', ...
            'rtfieldofview', 'rteffectivefocallength', 'rtpsf', ...
            'rtpsfdata', 'rtpsfwavelength', 'rtpsffieldheight', ...
            'rtpsfsamplespacing', 'rtpsfsupport', 'rtpsfsupportrow', ...
            'rtpsfsupportcol', 'rtotfdata', 'rtrelillum', ...
            'rtrifunction', 'rtriwavelength', 'rtrifieldheight', ...
            'rtgeometry', 'rtgeomfunction', 'rtgeomwavelength', ...
            'rtgeomfieldheight', 'rtgeommaxfieldheight'}          
        oType= 'optics';

    case {'chiefrayangle', 'chiefrayangledegrees', 'sensoretendue', ...
            'microlens', 'volts', 'digitalvalues', 'electrons', ...
            'dvorvolts', 'roielectrons', 'roivoltsmean', ...
            'roielectronsmean', 'hlinevolts', 'hlineelectrons', ...
            'vlinevolts', 'vlineelectrons', 'responseratio', ...
            'responsedr', 'analoggain', 'analogoffset', ...
            'sensordynamicrange', 'quantization', 'nbits', 'maxoutput', ...
            'quantizatonlut', 'quantizationmethod', ...
            'filtertransmissivities', 'infraredfilter', 'cfaname', ...
            'filternames', 'nfilters', 'filtercolorletters', ...
            'filtercolorletterscell', 'filterplotcolors', 'spectralqe', ...
            'pattern', 'dsnusigma', 'prnusigma', 'fpnparameters', ...
            'dsnuimage', 'prnuimage', 'columnfpn', 'columndsnu', ...
            'columnprnu', 'coloffsetfpnvector', 'colgainfpnvector', ...
            'noiseflag', 'reusenoise', 'noiseseed', 'pixel', ...
            'autoexpsoure', 'exposuretime', 'uniqueexptime', ...
            'exposureplane', 'cds', 'pixelvignetting', ...
            'nsamplesperpixel', 'sensormovement', 'movementpositions', ...
            'framesperpositions', 'sensorpositionsx', ...
            'sensorpositionsy', 'mccrecthandles', 'mcccornerpoints'}
        oType = 'sensor';
        
    case {'pdsize', 'fillfactor', 'pdarea', 'pdspectralqe', ...
            'conversiongain', 'voltageswing', 'wellcapacity', ...
            'darkcurrentdensity', 'darkcurrent', 'darkvoltage', ...
            'darkelectrons', 'readnoiseelectrons', 'readnoisevolts', ...
            'readnoisemillivolts', 'pdspectralsr', 'pixeldr'}
        oType = 'pixel';
        
    case {'render', 'colorbalance', 'colorbalancemethod', ...
            'demosaic', 'demosaicmethod', 'colorconversion', ...
            'colorconversionmethod', 'internalcolorspace', ...
            'internalcolormatchingfunciton', 'display', 'displayxyz', ...
            'displayxy', 'displaywhitepoint', 'displaymaxluminance', ...
            'displayspd', 'displaygamma', 'displaymaxrgb', ...
            'displaydpi', 'displayviewingdistance', 'l3'}
        oType = 'ip';
    otherwise
        % Maybe the default should be 'camera'?  Let's try.
        oType = 'camera';
        % error('Unparseable parameters %s\n', param);
        
end

end
