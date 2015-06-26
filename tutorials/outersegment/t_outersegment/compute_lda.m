function compute_lda(class1_L, class1_M, class1_S, class2_L, class2_M, class2_S)

% Still need to separate training and test data!

% all data
lda_data_all(:,1) = class1_L(:);
lda_data_all(:,2) = class1_M(:);
lda_data_all(:,3) = class1_S(:);

lda_data_all2(:,1) = class2_L(:);
lda_data_all2(:,2) = class2_M(:);
lda_data_all2(:,3) = class2_S(:);

tl = length(lda_data_all);
% training data
lda_datain(:,1) = lda_data_all(1:tl/2,1);
lda_datain(:,2) = lda_data_all(1:tl/2,2);
lda_datain(:,3) = lda_data_all(1:tl/2,3);

lda_datain2(:,1) = lda_data_all2(1:tl/2,1);
lda_datain2(:,2) = lda_data_all2(1:tl/2,2);
lda_datain2(:,3) = lda_data_all2(1:tl/2,3);

% test data
lda_data_test(:,1) = lda_data_all(tl/2+1:end,1);
lda_data_test(:,2) = lda_data_all(tl/2+1:end,2);
lda_data_test(:,3) = lda_data_all(tl/2+1:end,3);

lda_data_test2(:,1) = lda_data_all2(tl/2+1:end,1);
lda_data_test2(:,2) = lda_data_all2(tl/2+1:end,2);
lda_data_test2(:,3) = lda_data_all2(tl/2+1:end,3);

% figure; hold on;
% scatter3(cs1(:),cs2(:),cs3(:),'b');
% hold on; scatter3(csg1(:),csg2(:),csg3(:),'r');

lda_datain_all = [lda_datain; lda_datain2];
lda_type = [ones(1,length(lda_datain)) zeros(1,length(lda_datain2))];

lda_data_test_all = [lda_data_test; lda_data_test2];

linclass = ClassificationDiscriminant.fit(lda_datain_all,lda_type);

K = linclass.Coeffs(2,1).Const; % First retrieve the coefficients for the linear
L = linclass.Coeffs(2,1).Linear;% boundary between the second and third classes
                           % (versicolor and virginica).

% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2,x3) K + L(1)*x1 + L(2)*x2 +L(3)*x3;

for i = 1:length(lda_datain_all)
%     class_out(i) = f(lda_datain_all(i,1),lda_datain_all(i,2),lda_datain_all(i,3));
    class_out(i) = f(lda_data_test_all(i,1),lda_data_test_all(i,2),lda_data_test_all(i,3));
end

cl1end = length(lda_datain_all)/2;

[c1,c1x] = hist(class_out(1:cl1end),40);
[c2,c2x] = hist(class_out(cl1end+1:end),40);

% [c1,c1x] = hist(class_out(1:cl1end/2),40);
% [c2,c2x] = hist(class_out(cl1end+1:clend+1+(cl1end)/2),40);


% class_binary = class_out>-K;
class_unit_interval = (class_out)/(max(class_out)-min(class_out));
class_unit_interval = class_unit_interval-min(class_unit_interval);
% % class_unit_interval = class_unit_interval(class_unit_interval>-1.01);
% figure; hist(class_unit_interval)

[xm1 xm2] = meshgrid(20:2:50,20:2:50);
z = -(K+L(1)*xm1 + L(2)*xm2) ./L(3);
% figure; 
% hold on;
% surf(xm1,xm2,z);

% [tpr,fpr,thresholds] = roc(lda_type,class_out);

[tpr,fpr,thresholds] = roc(lda_type,class_unit_interval);
% figure; plotroc(lda_type,class_unit_interval);

[c1,c1x] = hist(class_unit_interval(1:cl1end),[0.05:0.05:1]);
[c2,c2x] = hist(class_unit_interval(cl1end+1:end),[0.025:0.05:1]);

figure;
subplot(121);
bar(c1x,c1./sum(c1),'barwidth',0.5);
hold on;
bar(c2x,c2./sum(c2),'r','barwidth',0.5);
my = max([c1./sum(c1) c2./sum(c2)]);
axis([0 1 0 1.1*my])
xlabel('Classifier output'); ylabel('Frequency');
% legend('Gabor present','Gabor absent');
% title(sprintf('CV = %d, Cont = %0.2f\nClassifier Output', cv, cont_mult));
title(sprintf('Classifier Output'));

subplot(122);
plot(fpr,tpr,'linewidth',2);
hold on;
plot(0.01:0.01:1,0.01:0.01:1,'k','linewidth',1);
% plotroc(lda_type,class_unit_interval);
fpr_diff = diff(fpr);
int_roc = sum(tpr(1:end-1).*fpr_diff);
% title(sprintf('CV = %d, Cont = %0.2f\nROC area = %0.2f',cv, cont_mult,int_roc));
title(sprintf('ROC area = %0.2f',int_roc));
xlabel('False Positive Rate'); ylabel('True Positive Rate');
toc

