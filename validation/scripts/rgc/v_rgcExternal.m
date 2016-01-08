clear

glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
movieFiles = dir([glmFitPath 'testmovie*']);
load([glmFitPath movieFiles(1).name], 'testmovie');
testmovieshort = testmovie.matrix(:,:,1:1+5*120);
% testmovieshort = (permute(testmovieshort, [2 1 3]));

scene = 0; sensor = 0; % not needed because using osIdentity object

os2 = osCreate('identity');
os2 = osSet(os2, 'rgbData', double(testmovieshort));

eyeAngle = 180; % degrees
eyeRadius = 3; % mm
eyeSide = 'right';

rgc2 = rgcPhys(scene, sensor, os2, eyeSide, eyeRadius, eyeAngle);
rgc2 = rgcSet(rgc2,'numberTrials',20);
rgc2 = rgcCompute(rgc2, os2);

rgc2linear = mosaicGet(rgc2.mosaic{1},'linearResponse');
rgc2psth = mosaicGet(rgc2.mosaic{1},'psthResponse');

% load('xvalall.mat');

client = RdtClient('isetbio');
client.credentialsDialog();
client.crp('resources/data/rgc')
[data, artifact] = client.readArtifact('xvalall', 'type', 'mat');
xvalall = data.xvalall;
%%
figure;
for i = 1%:36
    minlen = min([length(rgc2psth{i}) length(xvalall{i}.psth)]);
    diffpsth(i) = sum(abs(rgc2psth{i}(1:minlen) - xvalall{i}.psth(1:minlen)))./sum(.5*(rgc2psth{i}(1:minlen) + xvalall{i}.psth(1:minlen)));
    
%     subplot(6,7,i); hold on;
%     plot(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen),'b','linewidth',1);

    plot(rgc2psth{i}(1:minlen),'r ','linewidth',1);
    hold on;
    plot(xvalall{i}.psth(1:minlen),':k','linewidth',2);

    axis([0 6285 0  10]);%max(rgc2psth{i}(1:minlen))])
    [maxv, maxi] = max(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen)); title(sprintf('maxv = %.1f, maxi = %d',maxv,maxi));
    %     end;
end
% figure; plot(diffpsth,'x')

%%
% figure;
% for i = 1:39
%     subplot(6,7,i);
% %     plot(rgc2psth{i});
% %     hold on;
% %     plot(xvalall{i}.psth);
%     
%     plot(rgc2linear{i}(1:minlen/10)-xvalall{i}.lcif_const(1:10:minlen),'b','linewidth',2);
%     hold on;
% %     plot(xvalall{i}.cif0);
%     
% end
