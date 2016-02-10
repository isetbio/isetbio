


% figure; 
% hold on;
% nX = 0; nY = 0;
% % for cellTypeInd = 1:5
% cellTypeInd = 3;
% [nX,nY,~] = size(rgc1.mosaic{cellTypeInd}.linearResponse);
% % nX = nX + nXi;
% % nY = nY = nYi;
% t1 = (1:length(spPlot))';
% zfac = .0001;
% for xc = 1:nX
%     for yc = 1:nY
%         spPlot=rgc1.mosaic{cellTypeInd}.spikeResponse{xc,yc,1,2};
% %         spPlot=(median(horzcat(rgc1.mosaic{3}.spikeResponse{xc,yc,:,2})'));
%         plot3(zfac*(xc-ceil(nX/2))*t1+(xc-ceil(nX/2))+zeros(length(spPlot),1),1:length(spPlot),10*(yc-nY/2)+10*(yc-nY/2)*zfac*t1+spPlot)
%     end
% end
% view(0,0)
% view(-19,18);
% 
% %%
% for t = 1:10:6250
%     
% %     axis([-3 3 0+t 500+t -30 30])
%     axis([-4 4 0+t 1000+t -50 50])
%     pause(.005);
% end

%%

figure; 
% set(gcf,'position',[440   133   796   665]);
set(gcf,'position',[160 60 1070 740]);
hold on;
nX = 0; nY = 0;
% for cellTypeInd = 1:5
cellTypeInd = 3;
[nX,nY,~] = size(rgc1.mosaic{cellTypeInd}.linearResponse);
% nX = nX + nXi;
% nY = nY = nYi;
t1 = (1:length(spPlot))';
zfac = .002;
grid on;
for t = 1:100:5760
    hold on;
    cla
for xc = 1:nX
    for yc = 1:nY
        
        spPlot=rgc1.mosaic{cellTypeInd}.spikeResponse{xc,yc,1,2}(t:t+1000);
%         spPlot=(median(horzcat(rgc1.mosaic{3}.spikeResponse{xc,yc,:,2})'));

        t1 = (1:length(spPlot))';
        plot3(zfac*(xc-ceil(nX/2))*t1+(xc-ceil(nX/2))+zeros(length(spPlot),1),t/10000+1/10000*(1:length(spPlot)),10*(yc-nY/2)+10*(yc-nY/2)*zfac*t1+spPlot,'linewidth',2)
        
    end
end

xlabel('x position'); ylabel('time (sec)'); zlabel('y position');
set(gca,'fontsize',14);
% view(-19,18);

view(-2,2)

axis([-5 5 t/10000 (t+1000)/10000 -40 60]);
drawnow;
% pause(.01);
% hold off;
% clf

end
% view(0,0)
% view(-19,18);