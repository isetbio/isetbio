function [oType,p] = ieParameterOtype(param)
% Determine object type from parameter name
%
%   [oType,p] = ieParameterOtype(param,varargin)
%
% Parse the parameter string to determine which object (oType) and the
% parameter itself (p).
%
% Examples:
%   oType = ieParameterOtype('sensor dsnu')
%   [oType,p] = ieParameterOtype('optics/fnumber')
%
% Some will work without disambiguation because they are unique
%
%   oType = ieParameterOtype('dsnu')
%
% We also can take in the varargin, such as for units
%
%   [oType, p] = ieParameterOtype('pixel/size');
%
% This one won't work because size is available for too many objects
%
%    ieParameterOtype('size')
%
% This routine will require continual maintenace to maintain consistency
% with the specific object gets and sets.
%
% Examples:
%   [o p] = ieParameterOtype('dsnu sigma')
%   ieParameterOtype('sensor dsnu sigma')
%   [o p] = ieParameterOtype('oi size')
%   [o p] = ieParameterOtype('pixel size')
%   [o p] = ieParameterOtype('scene/hfov')
%   [o p] = ieParameterOtype('scene hfov')
%   [o p] = ieParameterOtype('optics/fnumber')
%
% Copyright Imageval LLC, 2013


%%
if ~exist('param','var') || isempty(param), error('Param required.'); end

%% If the param is simply one of the objects, return it as the oType

p = [];
switch ieParamFormat(param)
    case {'scene'}
        oType = 'scene'; return;
    case 'oi'
        oType = 'oi'; return;
    case 'optics'
        oType = 'optics'; return;
    case 'sensor'
        oType = 'sensor'; return;
    case 'pixel'
        oType = 'pixel'; return;
    case {'vci','ip'}
        oType = 'ip'; return;
end

%% Find the string before the first space or the first '/'.  

% Special characters are spaces and /
c1 = strfind(param,' ');   % Find the spaces
c2 = strfind(param,'/');   % Find the '/'

% Find the first space or '/'.  If we checked for whole string, we wouldn't
% need the bit above.  
pos = min([c1,c2]);

% Parse and return the string as oType
oType = [];
if ~isempty(pos)
    switch param(1:(pos-1))
        case 'scene'
            oType = 'scene';
        case 'oi'
            oType = 'oi';
        case 'optics'
            oType = 'optics';
        case 'sensor'
            oType = 'sensor';
        case 'pixel'
            oType = 'pixel';
        case {'vci','ip'}
            oType = 'ip';
    end
    
    % Check for success. Return the parameter, without the prepended term
    if ~isempty(oType), p = param((pos+1):end);  return; end
end

%% We didn't find the oType yet.

% Maybe param is one of the unique object parameter strings.  We have a
% list here.

% By the way ...
% I think there is a better way to do this in vistasoft, as per Adrian's
% coding using hashing.  Not sure we can use it in previous versions of
% Matlab, though.
p = param;
switch ieParamFormat(param)
    case {'objectdistance','meanluminance','luminance', ...
            'illuminant','illuminantname','illuminantenergy', ...
            'illuminantphotons','illuminantxyz','illuminantwave',...
            'illuminantcomment','illuminantformat'}
        oType = 'scene';
        
    case {'optics','opticsmodel','diffusermethod','diffuserblur'...
            'psfstruct','sampledrtpsf','psfsampleangles','psfanglestep',...
            'psfimageheights','raytraceopticsname','rtpsfsize',...
            }
        oType = 'oi';
        
    case {'fnumber','effectivefnumber','focallength','power',...
            'imagedistance','imageheight','imagewidth',...
            'numericalaperture','aperturedameter','apertureradius',...
            'magnification','pupilmagnification',...
            'offaxismethod','cos4thmethod','cos4thdata',...
            'otfdata','otfsize','otffx','otffy','otfsupport'...
            'psfdata','psfspacing','psfsupport',...
            'incoherentcutoffspatialfrequency','maxincoherentcutoffspatialfrequency',...
            'rtname','raytrace','rtopticsprogram','rtlensfile','rteffectivefnumber',...
            'rtfnumber','rtmagnification','rtreferencewavelength',...
            'rtobjectdistance','rtfieldofview','rteffectivefocallength',...
            'rtpsf','rtpsfdata','rtpsfsize','rtpsfwavelength',...
            'rtpsffieldheight','rtpsfsamplespacing',...
            'rtpsfsupport','rtpsfsupportrow','rtpsfsupportcol',...
            'rtotfdata','rtrelillum','rtrifunction','rtriwavelength',...
            'rtrifieldheight','rtgeometry','rtgeomfunction','rtgeomwavelength',...
            'rtgeomfieldheight','rtgeommaxfieldheight'}          
        oType= 'optics';
        
    case {'chiefrayangle','chiefrayangledegrees','sensoretendue',...
            'microlens','volts','digitalvalues','electrons',...
            'dvorvolts''roielectrons','roivoltsmean',...
            'roielectronsmean','hlinevolts','hlineelectrons',...
            'vlinevolts','vlineelectrons','responseratio','responsedr'...
            'analoggain','analogoffset','sensordynamicrange',...
            'quantization','nbits','maxoutput','quantizatonlut', ...
            'quantizationmethod','filtertransmissivities','infraredfilter',...
            'cfaname','filternames','nfilters','filtercolorletters',...
            'filtercolorletterscell','filterplotcolors','spectralqe',...
            'pattern','dsnusigma','prnusigma','fpnparameters',...
            'dsnuimage','prnuimage','columnfpn','columndsnu','columnprnu',...
            'coloffsetfpnvector','colgainfpnvector',...
            'noiseflag','reusenoise','noiseseed',...
            'pixel',...
            'autoexpsoure','exposuretime','uniqueexptime','exposureplane',...
            'cds','pixelvignetting',...
            'nsamplesperpixel'...
            'sensormovement','movementpositions','framesperpositions',...
            'sensorpositionsx','sensorpositionsy',...
            'mccrecthandles','mcccornerpoints'}
        oType = 'sensor';
        
    case {'pdsize','fillfactor','pdarea','pdspectralqe',...
            'conversiongain','voltageswing','wellcapacity',...
            'darkcurrentdensity','darkcurrent','darkvoltage','darkelectrons',...
            'readnoiseelectrons','readnoisevolts','readnoisemillivolts',...
            'pdspectralsr','pixeldr'}
        oType = 'pixel';
        
    case {'render','colorbalance','colorbalancemethod',...
            'demosaic','demosaicmethod','colorconversion','colorconversionmethod',...
            'internalcolorspace','internalcolormatchingfunciton',...
            'display','displayxyz','displayxy','displaywhitepoint',...
            'displaymaxluminance','displayspd','displaygamma','displaymaxrgb',...
            'displaydpi','displayviewingdistance','l3'}
        oType = 'ip';
    otherwise
        % Maybe the default should be 'camera'?  Let's try.
        oType = 'camera';
        % error('Unparseable parameters %s\n',param);
        
end

end

