function mosaicAverageGLM = mosaicAverage(mosaicGLM)

for i = 1:length(mosaicGLM)
    
    sRFtemp = mosaicGLM{i}.linearfilters.Stimulus.space_rk1;
    sRF(i,:,:) = sRFtemp;
%     sRF(i,:) = reshape(sRFtemp, size(sRFtemp,1)*size(sRFtemp,2));
    maxsRF(i)  = max(sRFtemp(:));
    minsRF(i)  = min(sRFtemp(:));
    meansRF(i) = mean(sRFtemp(:));
    
    tCtemp = mosaicGLM{i}.linearfilters.Stimulus.time_rk1;
    tC(i,:) = tCtemp;
    maxtC(i) = max(tCtemp);
    mintC(i) = min(tCtemp);
    meantC(i) = mean(tCtemp);
    
    tonicD(i,:) = 0;% mosaicGLM{i}.linearfilters.tonicDrive{i,1}(:);
    
    nlcoeffs(i,:) = mosaicGLM{i}.model.Coefficients.Estimate;
    if isfield(mosaicGLM{i}.stafit,'center_sd_x')
    sd_x(i,:) = mosaicGLM{i}.stafit.center_sd_x; 
    sd_y(i,:) = mosaicGLM{i}.stafit.center_sd_y;
    end
end

% meanRF = mean(sRF);
% sRFrs = reshape(meanRF,size(sRFtemp,1),size(sRFtemp,2));
% figure; surf(sRFrs);

oldSize = size(squeeze(sRF(i,:,:)),1);
newSize = size(sRFtemp,1)*3+1;
meanRF = zeros(newSize);
meanCtr = zeros(newSize);

for i = 1:length(mosaicGLM)
    [maxPR(i) maxPRind(i)] = max(max(abs(squeeze(sRF(i,:,:))),[],1));
    [maxPC(i) maxPCind(i)] = max(max(abs(squeeze(sRF(i,:,:))),[],2));
   
    xv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - (maxPCind(i) - floor(oldSize/2)) ;
    yv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - (maxPRind(i) - floor(oldSize/2));
    
    meanRF(xv,yv) = meanRF(xv,yv) + squeeze(sRF(i,:,:));
    meanCtr(xv,yv) = meanCtr(xv,yv) + ones(size(squeeze(sRF(i,:,:))));
end

% figure; surf(meanRF); shading flat

meanAvg = meanRF./meanCtr;
meanAvg(meanCtr<4) = 0;
% figure; surf(meanAvg); shading flat; shading interp

sum(meanAvg(:));%  3.8475 
sum(sRF(:))/length(mosaicGLM); % 4.5199

max(meanAvg(:)); % 0.2039
mean((max(reshape(sRF,size(sRF,1), size(sRF,2)*size(sRF,3))'))) ;% 0.2738

min(meanAvg(:)); %  -0.0377
mean((min(reshape(sRF,size(sRF,1), size(sRF,2)*size(sRF,3))'))); % -0.0459


% figure; surf(meanCtr);
% 
% figure; plot(mean(tC));

max(mean(tC)); % 1.119floor(oldSize/2)
mean(max(tC'));% 1.1279

mean(tC(:)); % 0.0243

min(mean(tC));% -0.2837

mean(min(tC')); %  -0.29floor(oldSize/2)3

mean(tonicD);% 2.2702

cv = [(floor(newSize/2)+1) - floor(oldSize/2) : (floor(newSize/2)+1) + floor(oldSize/2)] - 1;
mosaicAverageGLM.linearfilters.Stimulus.space_rk1 = meanAvg(cv,cv);
mosaicAverageGLM.linearfilters.Stimulus.time_rk1 = mean(tC);
mosaicAverageGLM.modelavg = mean(nlcoeffs);

% Need to change this
mosaicAverageGLM.model = mosaicGLM{1}.model;

mosaicAverageGLM.sd = [sqrt(mean(sd_x(find(sd_x~=0)).^2)) sqrt(mean(sd_y(find(sd_y~=0)).^2))];
