function mosaicPlot(ir, bp, sensor, irparams, cellType, ecc)
% Plots the cone, bipolar and RGC mosaic on top of one another
% 
%       mosaicPlot(ir, bp, sensor, irparams, cellType, ecc)
% 
% 
% (c) 2016 JRG isetbio team

figure;

spacing = irGet(ir,'spacing');

% 
cone_mosaic = sensor.human.coneType;

scale2 = size(cone_mosaic,2)/spacing;%1;%size(cone_mosaic,2)/40;
scale1 = scale2;
  
% [xg yg] = meshgrid([1:size(cone_mosaic,1),1:size(cone_mosaic,1)]);
% xg2 = xg(1:size(cone_mosaic,1),1:size(cone_mosaic,2)); yg2 = yg(1:size(cone_mosaic,1),1:size(cone_mosaic,2));
% 
% % figure; scatter(xg2(:),yg2(:),40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])
% % figure; 
% hold on;
% scatter(yg2(:)./scale1,xg2(:)./scale2,40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])

% figure; 
% imagesc(4-cone_mosaic); colormap jet;
imagesc((1:size(cone_mosaic,2))./scale1, (1:size(cone_mosaic,1))./scale1, 4-cone_mosaic); colormap jet;

    % circle samples
bpRad = 1;%size(bp.sRFcenter,1);
circle_samples = 0:0.05:2*pi;
x_circle = bpRad*cos(circle_samples);
y_circle = bpRad*sin(circle_samples);
% 
szbp = size(bp.sRFcenter,1)-.5;
% figure;
hold on;
szmosaic = size(cone_mosaic);
ycind = 0; xcind = 0;
for yc = 1:szbp*2:szmosaic(1)
    ycind = ycind+1;
    for xc = 1:szbp*2:szmosaic(2)
        xcind = xcind+1;
%         plot((szbp*y_circle+yc)./scale2,(szbp*x_circle+xc+(szbp)*mod(ycind+1,2))./scale1,'linewidth',1,'color','c');
        plot((szbp*x_circle+xc+0*(szbp)*mod(ycind+1,2))./scale1,(szbp*y_circle+yc)./scale2,'linewidth',1,'color','c');
    end
end


%
for i = 1:length(ir.mosaic{1}.cellLocation); 
    loc(i,:) = irparams.inputScale.*ir.mosaic{1}.cellLocation{i}; 
end;
% figure; scatter(loc(:,1),loc(:,2),100); axis equal

rfRad = 1*ir.mosaic{1}.rfDiaMagnitude/2;
circle_samples = 0:0.05:2*pi;
c1(1,:) = rfRad*cos(circle_samples);
c1(2,:) = rfRad*sin(circle_samples);

% figure; 
hold on;
for i = 1:length(ir.mosaic{1}.cellLocation)
    plot((loc(i,2)+c1(1,2:end)-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(1))./scale1,(loc(i,1)+c1(2,2:end)-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(2))./scale2,'m','linewidth',3); axis equal
%         plot(loc(i,2)+c1(1,2:end)+1*rfRad-1*ir.mosaic{1}.cellLocation{1,1}(1),loc(i,1)+c1(2,2:end)+rfRad-2*ir.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,1)+c1(1,2:end)+2*rfRad-1*ir.mosaic{1}.cellLocation{1,1}(1),loc(i,2)+c1(2,2:end)+rfRad-1*ir.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,2)+c1(2,2:end)-16.5,loc(i,1)+c1(1,2:end)-16.5,'m','linewidth',3); axis equal
    
end
axis equal
% axis([1 94 1 74])
title(sprintf('%s simulated mosaic at %1.1f\\circ Ecc',cellType(1:end-3),ecc/.3));
set(gca,'fontsize',16)
xlabel(sprintf('Distance (\\mum)')); 
ylabel(sprintf('Distance (\\mum)'));