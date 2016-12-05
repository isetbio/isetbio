function timeAxis = rgcInitTime(rgcM, innerRetina, mosaicInd)
% Initialize a temporal impulse response function for a mosaic
%
%   timeAxis = rgcInitTime(rgcM, innerRetina, mosaicInd)
% 
% Build a temporal impulse response function for an rgc mosaic
%
% Inputs:
%   rgcM:         The calling RGC Mosaic object
%   innerRetina:  The parent inner retina object 
%   mosaicInd:    Which mosaic in innerRetina.mosaic{} we are initializing
%
% Temporal impulse response peak magnitude scale parameters. see "Spatial
% Properties and Functional Organization of Small Bistratified Ganglion
% Cells in Primate Retina", Field, et al., J. Neuroscience, 2007, Fig. 1.
%
% See also:
%  
%
% JRG/BW ISETBIO Team, 2016

%%  Do we allow other parameters?
p = inputParser;
p.addRequired('rgcM')
p.addRequired('innerRetina')
p.addRequired('mosaicInd',1,@isscalar);

%%
rfTempMult = [1 -1 1 -1 1];       % invert sign for OFF types

% The temopral impulse response function for each cell consists of three
% vectors, one for each RGB channel. They have similar shapes but different
% polarities and relative magnitudes depending on the cell type.
switch ieParamFormat(rgcM.cellType)
    case{'smallbistratified','sbc'}
        % SBCs have the B channel impulse response reversed in polarity
        rgbTempMult = [-0.4 -0.4 1]; % [R G B] magnitudes
    otherwise
        % The other four of the big five types have the same polarity
        rgbTempMult = [0.4 1 0.4]; % [R G B] magnitudes
end

% Get the sampling interval set in the outersegment or sensor
integrationTime = innerRetina.timing;
for rgbInd = 1:3
    % scale for differences in RGB channel and ON/OFF type
    multFactor = rgbTempMult(rgbInd)*rfTempMult(mosaicInd);
    % Build impulse responses for center and surround; usually the same.
    params.filterDuration = 0.2;  % 200 ms
    params.samplingTime   = innerRetina.timing;  % Time step.  Usually 1 or 2 ms
    [tmp, timeAxis] = rgcImpulseResponsePillow(params);
    % Old code:
    % [tmp, timeAxis] = buildTemporalImpulseResponse(integrationTime);
    
    rgcM.tCenter{rgbInd,1} = multFactor*tmp;
    rgcM.tSurround{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    
    % vcNewGraphWin; plot(timeAxis,tmp);
end

end