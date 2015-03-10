function ISA = sensorVignetting(ISA,pvFlag,nAngles)
%Compute etendue (vignetting) slot and place in sensor structure
%
%   ISA = sensorVignetting(ISA,pvFlag)
%
%Examples:
%   foo = sensorVignetting; plotSensorEtendue(foo);
%   foo = sensorVignetting([],3); plotSensorEtendue(foo);
%
%   ISA = vcGetObject('ISA'); ISA = sensorSet(ISA,'pixelVignetting',1);
%   foo = sensorVignetting(ISA); plotSensorEtendue(foo);
%
% Copyright ImagEval Consultants, LLC, 2006.

if notDefined('ISA'),    ISA = vcGetObject('ISA'); end
if notDefined('pvFlag'), pvFlag = sensorGet(ISA,'pixelVignetting'); end
if notDefined('nAngles'),nAngles=5; end

if isempty(pvFlag), pvFlag = 0; end

switch pvFlag
    case 0   %skip, so set etendue to 1s
        sz = sensorGet(ISA,'size');
        ISA = sensorSet(ISA,'etendue',ones(sz));
    case 1   %bare, nomicrolens
        ISA = mlAnalyzeArrayEtendue(ISA,'nomicrolens',nAngles);
    case 2   %centered
        ISA = mlAnalyzeArrayEtendue(ISA,'centered',nAngles);
    case 3  %optimal
        ISA = mlAnalyzeArrayEtendue(ISA,'optimal',nAngles);
    otherwise
        error('Unknown pvFlag %s\n',pvFlag);
end

end