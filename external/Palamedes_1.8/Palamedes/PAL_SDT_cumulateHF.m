%
%PAL_SDT_cumulateHF takes a matrix of hits and false alarms from a 
%rating scale experiment and returns matrices with the cumulative number 
%of hits and false alarms, the total number of signal and noise trials 
%(or signal 1 and signal 2 trials), and the probabilities of hits and false 
%alarms. The cumulative number of hits and false alarms for the last pair 
%of entries, which equals the total numbers of trials, is ommitted from
%the output.
%
%The routine can be used to provide the input arguments for
%PAL_SDT_ROCML_Fit
%
%Internal function
%
%Introduced: Palamedes version 1.6.0 (FK & NP)

function [cumNumHF, OutOfNum, pHF] = PAL_SDT_cumulateHF(NumHF)

NumTrials = sum(NumHF);
cumNumHF = cumsum(NumHF);
cumNumHF = cumNumHF(1:size(cumNumHF)-1,:);
OutOfNum = repmat(NumTrials,length(cumNumHF),1);
pHF = cumNumHF./OutOfNum;
