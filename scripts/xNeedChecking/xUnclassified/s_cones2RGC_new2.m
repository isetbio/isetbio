
clear

load('sensor_new_CV_2_Cont_0020.mat')

noiseflag = 1;
linearOS = osLinear('noiseFlag',noiseflag);
% adaptedOS = osBioPhys('noiseFlag',noiseflag);
% compute outersegment current output
linearOS = osLinearCompute(linearOS,sensor,params);

%% Create RGC structure

% creating the RGC parameter object
rgcP = rgcParameters;

% The absorptions structure is also the RGC parameter data
rgcP.set('cone voltages',linearOS.ConeCurrentSignal); % this should be cone current
rgcP.set('sensor',sensor);
rgcP.set('oi',oi);

% What is the default?  RGC spacing the same as cone spacing (I think).
rgcP.addLayer('on parasol', 20);  
rgcP.addLayer('off parasol');  
 
rgcP.addLayer('on midget');  
rgcP.addLayer('off midget');  

rgcP.addLayer('small bistratified');  % Sign of signal seems wrong. Check.

%% RGC Spike computation


nL = rgcP.get('nL');
layer = cell(1,nL);
for ii=1:nL, layer{ii} = rgcP.get('layer',ii); end

hasFeedback = ones(1,nL);
for ii=1:nL, layer{ii}.hasFeedback = hasFeedback(ii); end

hasCoupling = zeros(1,nL);
for ii=1:nL, layer{ii}.hasCoupling = hasCoupling(ii); end

% Spike threshold as a percentage of the RGC voltage swing
vSwingOriginal = layer{1}.get('vswing');
for ii=1:nL, layer{ii}.set('rgc volt thresh',0.2*vSwingOriginal); end

% layer{1}.set('cone weights',[.3 0 0; 0 0 .3]);
% layer{1}.get('cone weights')

% We also want the largest value in the coupling and feedback to be smaller
% than the spike threshold.  So attend to that around here.
% fbTimeOriginal = layer{1}.get('fbtr');
% for ii=1:nL
%     layer{ii}.set('rgc volt thresh',0.5*fbTimeOriginal);
% end

rgcComputeSpikes(rgcP);


%%
figure; 
for spind = 1:4
    subplot(2,2,spind); 
    imagesc(rgcP.layers{1,spind}.rgcvTimeSeries);
    
    switch spind
        case 1
            title('on parasol');
        case 2
            title('off parsol');
        case 3
            title('on midget');
        case 4 
            title('off midget');
    end
xlabel('time (ms)'); ylabel('RGC Spikes');
end
