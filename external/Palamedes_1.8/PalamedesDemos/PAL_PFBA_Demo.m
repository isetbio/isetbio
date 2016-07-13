%
%PAL_PFBA_Demo  Demonstrates use of Palamedes routine PAL_PFBA_Fit to fit a
%Psychometric Function (PF) to some data using a Bayesian criterion and
%determine the standard error estimates of the PF's parameters.
%
%Demonstrates usage of Palamedes function:
%-PAL_PFBA_Fit
%secondary:
%PAL_Logistic
%PAL_pdfNormal
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFBA_Fit
%
%NP (September 2009)

clear all;  %Clear all existing variables from memory

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure does not render exactly as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

%Stimulus intensities
StimLevels = [-3:1:3];

%Number of positive responses (e.g., 'yes' or 'correct' at each of the 
%   entries of 'StimLevels'  
NumPos = [55 55 66 75 91 94 97];

%Number of trials at each entry of 'StimLevels'
OutOfNum = [100 100 100 100 100 100 100];

%Define parameter space:
priorAlphaValues = [-1:.01:1];
priorBetaValues = [-.5:.01:.5];

%Set Guess-rate and Lapse-rate (fixed parameters)
gamma = 0.5;
lambda = 0;

%Fit a logistic function
PF = @PAL_Logistic;

%Define a prior distribution across parameter space (optional)
[Alpha Beta] = meshgrid(priorAlphaValues,priorBetaValues);
prior = PAL_pdfNormal(Alpha,0,1).*PAL_pdfNormal(Beta,0,1);
prior = prior./sum(sum(prior));

%Fit function
[paramsValues posterior] = PAL_PFBA_Fit(StimLevels, NumPos, OutOfNum, ...
    priorAlphaValues, priorBetaValues, gamma, lambda, PF, 'prior',prior);

%Put summary of results on screen
message = sprintf('\rThreshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f',paramsValues(2));
disp(message);
message = sprintf('Standard error of Threshold: %6.4f',paramsValues(3));
disp(message);
message = sprintf('Standard error of Slope: %6.4f',paramsValues(4));
disp(message);

%Create simple contour plot of posterior
figure('name','Bayesian Psychometric Function Fitting');
contour(posterior)
set(gca, 'Xtick',[1:50:201],'FontSize',12);
set(gca, 'XtickLabel', {'-1','-.5','0','.5','1'});
xlabel('Alpha');
set(gca, 'Ytick',[1:25:101],'FontSize',12);
set(gca, 'YtickLabel', {'-.5','-.25','0','.25','.5'});
ylabel('Log(Beta)');