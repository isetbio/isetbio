% Illustrate noise generation for the coneMosaic object
%
% Description:
%    Illustrates the various ways to generate noise using the coneMosaic
%    object. We set coneMosaic.noiseFlag to 'random', 'frozen', and 'none'.
%
%    This tutorial also illustrates how to access some of the key fields of
%    the coneMosaic object programatically and make some plots of what is
%    in there.
%
%    It does not yet illustrate controlling the 'seed' for the 'frozen'
%    noise case.
%
% See Also:
%   t_conesMosaicBasic.
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    08/04/17  dhb  Clean up.
%    07/25/18  jnm  Formatting

%% Initialize
ieInit;

%% Create a simple scene and compute the optical image
% In ISETBio, Poisson noise isn't added to the retinal image by default,
% instead we take the mean value and allow the noise emerge whene we
% compute isomerizations.
scene = sceneCreate('slanted bar');
oi = oiCreate;
oi = oiCompute(oi, scene);

%% Look at effect of random Poisson isomerization noise.
% Random meaning not frozen -> different on every run.
%
% Create a cone mosaic and generate isomerizations. The noiseFlag is
% 'random' by default, but it is set here explictily just because.
cMosaic = coneMosaic;
cMosaic.noiseFlag = 'random';
cMosaic.compute(oi);

%% View the mosaic in a window. Note that the isomerizations look noisy.
% If you want, select Plot:hLine LMS in the window and click on a row. Up
% will come a plot showing the L, M and S cone isomerizations in separate
% subplots, and you will observe that they look noisy. Just below, we'll
% make that type of plot programatically.
cMosaic.window;

%% Make a line plot of L, M, and S isomerizations.
% Here we will show how to make the plot without using the window, but
% instead reaching down into the object to get what we need and then
% makeing the plot directly.
%
% We may eventually incorporate this functionality into coneMosaic.plot,
% but it is not there currently. The code below was adopted from the
% plotting engine used by coneMosaic.window, look in file coneMosaicWindow.
%
% The pattern property of the coneMosaic tells us which type of cone is at
% each position, where 2 -> L, 3 -> M, 4 -> S. (Locations with 1 have no
% cone). The code loops over each cone type, pulls out the isomerizatoins
% for each cone type, and plots each in a separate subpanel.
%
% Note that the field that contains the isomerizations is called
% 'absorptions'. In ISETBio, absorptions is a synonym for isomerizations.
linePlot1 = vcNewGraphWin([], 'tall');
whichLine = 15;
yLabelConeNames = ['L', 'M', 'S'];
yLabelCommon = 'Cone Isomerizations';
plotSpec = {'ro-', 'go-', 'bo-'};
for ii = 2:4
    subplot(3, 1, ii - 1);
    positionIndex = find(cMosaic.pattern(whichLine, :) == ii);
    plot(positionIndex, cMosaic.absorptions(whichLine, positionIndex), ...
        plotSpec{ii - 1}, 'LineWidth', 2);
    grid on;
    xlabel('Horizontal Position (cones)');
    ylabel([yLabelConeNames(ii - 1) ' ' yLabelCommon]);
    set(gca, 'xlim', [1 size(cMosaic.absorptions, 2)]);
end

%% Do a second draw and add to the plot
% This lets you see that the noise is different on the second draw.
%
% Note that the S cone isomerizations appear noisier. If you look at the
% y-axis scale in the three plots, you'll see that there are fewer S cone
% isomerizations. This is primarily because of wavelength-dependent
% absorption by the lens and macular pigments. The smaller number of
% isomerizations means that the Poisson noise has a larger relative effect
% for the S cones.
cMosaic.compute(oi);
figure(linePlot1);
plotSpec = {'ro:', 'go:', 'bo:'};
for ii = 2:4
    subplot(3, 1, ii - 1);
    hold on
    positionIndex = find(cMosaic.pattern(whichLine, :) == ii);
    plot(positionIndex, cMosaic.absorptions(whichLine, positionIndex), ...
        plotSpec{ii - 1}, 'LineWidth', 2);
end
title('Noisy Isomerizations (Two Separate Draws)');

%% Generate isomerizations with no noise
% This returns the mean number of isomerizations, without Poisson noise.
% Look how beautifully smooth these functions are.
cMosaic.noiseFlag = 'none';
cMosaic.compute(oi);
linePlot2 = vcNewGraphWin([], 'tall');
whichLine = 15;
yLabelConeNames = ['L', 'M', 'S'];
yLabelCommon = 'Cone Isomerizations';
plotSpec = {'ro-', 'go-', 'bo-'};
for ii = 2:4
    subplot(3, 1, ii - 1);
    positionIndex = find(cMosaic.pattern(whichLine, :) == ii);
    plot(positionIndex, cMosaic.absorptions(whichLine, positionIndex), ...
        plotSpec{ii - 1}, 'LineWidth', 2);
    grid on;
    xlabel('Horizontal Position (cones)');
    ylabel([yLabelConeNames(ii-1) ' ' yLabelCommon]);
    set(gca, 'xlim', [1 size(cMosaic.absorptions, 2)]);
end
title('Noise Free (Mean) Isomerizations (Two Separate Draws)');

% Save the mean isomerizations from the no noise case above for plotting.
meanIsomerizations = cMosaic.absorptions;

%% Frozen noise
% The isomerizations are noisy, but the same every time. That's shown in
% the plot where the first and second draws are plotted against the
% corresponding mean isomerization for all cones. Each red dot (first
% draw) has a black dot (second draw) right on top of it.
%
% Also note how variance increases with isomerizations, as it should for
% Poisson noise.
%
cMosaic.noiseFlag = 'frozen';
cMosaic.compute(oi);
frozenIsomerizations1 = cMosaic.absorptions;
cMosaic.compute(oi);
frozenIsomerizations2 = cMosaic.absorptions;
vcNewGraphWin;
hold on
plot(meanIsomerizations(:), frozenIsomerizations1(:), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(meanIsomerizations(:), frozenIsomerizations2(:), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 4);
xlabel('Mean Isomerizations');
ylabel('Frozen Isomerizations');
title('Two Draws of Frozen Isomerizations Against Mean');
legend({'First Frozen Draw', 'Second Frozen Draw'}, ...
    'Location', 'NorthWest');
axis('square');
