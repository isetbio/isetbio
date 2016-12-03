function timeAxis = rgcInitTime(rgcM, innerRetina, mosaicInd)
% Generate temporal impulse response functions for R, G, B channels
%
%
% 

%%
% ir = rgcM.parent;


%%  Do we allow other parameters?

% Temporal impulse response peak magnitude scale parameters.
% see "Spatial Properties and Functional Organization of Small
% Bistratified Ganglion Cells in Primate Retina", Field, et al.,
% J. Neuroscience, 2007, Fig. 1.

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
    % Build the separate impulse responses for center and surround; usually
    % the same.
    [tmp, timeAxis] = buildTemporalImpulseResponse(integrationTime);
    rgcM.tCenter{rgbInd,1} = multFactor*tmp;
    rgcM.tSurround{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    % plot(timeAxis,tmp);
end

end