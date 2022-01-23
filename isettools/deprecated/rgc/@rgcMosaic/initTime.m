function timeAxis = initTime(rgcM, innerRetina)
% Initialize a temporal impulse response function for a mosaic
%
% Syntax:
%   timeAxis = rgcMosaic.initTime(rgcM, innerRetina)
%
% Description:
%    Build a temporal impulse response function for an rgc mosaic.
%
% Inputs:
%    rgcM        - Object. The calling RGC Mosaic object
%    innerRetina - Object. The parent inner retina object
%    mosaicInd   - Numeric. The index of which mosaic in
%                  innerRetina.mosaic{} to initialize.
%
% Outputs:
%    timeAxis    - Array. The time axis for the temporal impulse response.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * Temporal impulse response peak magnitude scale parameters. see
%      "Spatial Properties and Functional Organization of Small
%      Bistratified Ganglion Cells in Primate Retina", Field, et al., J.
%      Neuroscience, 2007, Fig. 1.
%

% History:
%    XX/XX/16  JRG/BW ISETBIO Team, 2016
%    06/14/19  JNM    Documentation pass

%%  Do we allow other parameters?
p = inputParser;
p.addRequired('rgcM')
p.addRequired('innerRetina')
p.parse(rgcM, innerRetina);

%% Use the Pillow function, I guess.
% Build impulse responses for center and surround; usually the same.
% But maybe we should have different impulse for center and surround?
params.filterDuration = 0.2;  % 200 ms
params.samplingTime = innerRetina.timeStep;  % Time step. Usually 1 or 2 ms
[tmp, timeAxis] = rgcImpulseResponsePillow(params);

rgcM.tCenter = tmp;
rgcM.tSurround = -1 * tmp;
% vcNewGraphWin;
% plot(timeAxis, tmp);

% Old code:
% [tmp, timeAxis] = buildTemporalImpulseResponse(integrationTime);

end

% We had this complicated code for creating the temporal impulse responses,
% one for each RGB (I think LMS was meant). And then there were scale
% factors depending on RGB. This may have had something to do with the
% Chichilnisky lab analysis that JRG put in here.
%
% The temopral impulse response function for each cell consists of three
% vectors, one for each RGB channel. They have similar shapes but different
% polarities and relative magnitudes depending on the cell type.
%
% switch ieParamFormat(rgcM.cellType)
%     case{'smallbistratified', 'sbc'}
%         % SBCs have the B channel impulse response reversed in polarity
%         rgbTempMult = [-0.4 -0.4 1]; % [R G B] magnitudes
%     otherwise
%         % The other four of the big five types have the same polarity
%         rgbTempMult = [0.4 1 0.4]; % [R G B] magnitudes
% end
