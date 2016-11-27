% t_osLinearize
%
%   Illustrates the foveal and peripheral cone outer segment temporal
%   impulse responses using the osBioPhys object at different mean
%   intensities.
%
% JRG/BW 11/2016 (c) isetbio team

%%
ieInit;

cMosaic = coneMosaic('os',osLinear,'pattern',[2 2 2]);

%% Background levels

% Mean of isomerization stimulus in R*/sec (scaled by time step to be
% placed in bins of photons). The mean rate affects the magnitude of
% the impulse response due to adaptation.
meanIsoArray = [100 500 1000 2000 5000 10000]*cMosaic.integrationTime;

%%  Loop on different background rates and plot
fovea = zeros(4002,length(meanIsoArray));
periphery = zeros(4002,length(meanIsoArray));

for ii = 1:length(meanIsoArray)
    
    % Set mean background for start of current output
    cMosaic.absorptions = repmat(meanIsoArray(ii),1,3);
        
    % Compute outer segment currents.
    tmp     = os.linearFilters(cMosaic,'eccentricity',true);
    fovea(:,ii) = tmp(:,1);
    
    tmp = os.linearFilters(cMosaic,'eccentricity',false);
    periphery(:,ii) = tmp(:,1);
    
end

%%
vcNewGraphWin([],'wide'); 
T = os.timeAxis;
subplot(1,2,1); hold on;
plot(T,fovea,'LineWidth',3);
grid on; % legend('Peripheral','Foveal');

% axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
xlabel('Time (sec)','FontSize',14);
ylabel('pA','FontSize',14);
set(gca,'fontsize',14);
title(sprintf('Foveal IR (R*/sec)'));

%
subplot(1,2,2); hold on;
plot(T,periphery,'LineWidth',3);
grid on; % legend('Peripheral','Foveal');
% axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
xlabel('Time (sec)','FontSize',14);
ylabel('pA','FontSize',14);
set(gca,'fontsize',14);
title(sprintf('Peripheral IR (R*/sec)'));

% legend(num2str(meanIsoArray));
legend('100', '500', '1000', '2000', '5000', '10000')
