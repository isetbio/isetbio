%
%PAL_SDT_Summ_MultiplePF_SimulateObserverParametric  Simulate observer
%for a summation experiment with multiple PFs (psychometric functions)
% based on SDT (signal detection theory)
%
%Syntax: 
%
% NumPos = PAL_SDT_Summ_MultiplePF_SimulateObserverParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q)
%
%Input:
%
%   'StimLevels': matrix of stimulus levels
%       
%   'gParams' vector for stimulus level scaling factors gA, gB gC etc. 
%       These would normally be the params fitted to the PF data
%
%   'pParams' vector for corresponding transducer exponents pA, pB, pC etc.
%       These would normally be the params fitted to the PF data
%   
%   'OutOfNum': number of trials per stimulus level
%
%   'SummFunc': either @PAL_SDT_PS_uneqSLtoPC, @PAL_SDT_AS_uneqSLtoPC
%                       @PAL_SDT_PS_SLtoPC, @PAL_SDT_AS_SLtoPC;
%
%   'M': number of alternatives/intervals in forced-choice task
%
%   'Q': number of monitored channels
%
%Output:
%
%   'NumPos'  simulated number of correct reponses for each stimulus level
%
%
%Example:
%
% This example uses the routine to generate simulated numbers of correct
% detections based on a probability summation model fitted simultaneously 
% to three psychometric functions, one for stimulus A, one for stimulus B 
% and one for stimulus A+B, under the Fixed Attention Window scenario, i.e. 
% when the observer is monitoring both A and B channels.  The summation model 
% allows for unequal stimulus levels for A and B 
%
% StimLevels(1,1,:)=[1 2 3 4 5 6 7 8]; % for stim A alone
% StimLevels(1,2,:)=[0 0 0 0 0 0 0 0]; % channel B for stim A alone
% StimLevels(2,1,:)=[0 0 0 0 0 0 0 0]; % channel A for stim B alone
% StimLevels(2,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B alone
% StimLevels(3,1,:)=[1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2]; %Stim A in Stim A+B
% StimLevels(3,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B in Stim A+B
% 
% OutOfNum(1,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(2,:) = [80 80 80 80 80 80 80 80];
% OutOfNum(3,:) = [80 80 80 80 80 80 80 80];
% 
% SummFunc=@PAL_SDT_PS_uneqSLtoPC; probability summation (PS) model for
%       unequal stimulus levels for A and B
%
% M=2;  % 2-AFC
% Q=2;  % two monitored channels
%
% gParams=[0.2998 0.2508]; %g parameters for PS model
% pParams=[1.2499 1.4779]; %p parameters for PS model
%
%Obtain simulated number positive responses:
%
%   NumPos = PAL_SDT_Summ_MultiplePF_SimulateObserverParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q)
%
%returns something like:
%
% NumPos =
%
%     45    47    49    62    59    69    72    77
%     41    50    51    59    68    69    77    76
%     46    48    55    69    74    77    78    81
%
% 
%Introduced: Palamedes version 1.8.0 (FK & NP)

function NumPos = PAL_SDT_Summ_MultiplePF_SimulateObserverParametric(StimLevels,gParams,pParams,OutOfNum,SummFunc,M,Q)

[nrows, ncols, nnums]=size(StimLevels);

PC=zeros(nrows,nnums);
NumPos=zeros(nrows,nnums);
StimVec=zeros(1,ncols); % vector for corresponding stimulus levels in summation equation

for j=1:nrows
    for i=1:nnums
        for k=1:ncols
            StimVec(1,k)=StimLevels(j,k,i);
        end
        PC(j,i) = SummFunc(StimVec,gParams,pParams,M,Q);
    end
end


for j=1:nrows
    for i=1:nnums
        Pos = rand(OutOfNum(j,i),1);
        Pos(Pos < PC(j,i)) =1;
        Pos(Pos ~= 1) = 0;
        NumPos(j,i) = sum(Pos);
    end
end

