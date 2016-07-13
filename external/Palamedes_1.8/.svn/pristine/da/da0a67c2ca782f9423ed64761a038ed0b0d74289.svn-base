%
%PAL_SDT_PF_SimulateObserverParametric  Simulate observer characterized by
%   a PF (psychomtric function) from an SDT (signal detection thory) model.
%
% NumPos = PAL_SDT_PF_SimulateObserverParametric(StimLevels,params,...
%   OutOfNum,SDTfunc,M)
%
%Input:
%
%   'StimLevels': vector of stimulus levels
%   
%   'params': 2-valued vector of parameters. These are the stimulus level 
%       scaling factor g and the transducer exponent p.
%
%   'OutOfNum': number of trials per stimulus level
%
%   'SDTfunc': examples @PAL_SDT_2AFC_DPtoPC, 
%       @PAL_SDT_3AFCoddity_IndMod_DPtoPC;
%
%   'M': number of alternatives in forcd-choice task.  If SDTfunc does not
%       take M as an argument M should be empty, i.e. 'M=[]'
%
%Output:
%
%   'NumPos': vector of number of simulated postive responses
%
%Example:
%
%   StimLevels=[1 2 3 4 5 6];  
%   params=[0.3222 1.8320];
%   OutOfNum=[100 100 100 100 100 100];
%   SDTfunc=@PAL_SDT_2AFC_DPtoPC;
%   M=[];
%
%Simulate data:
%
%   NumPos = PAL_SDT_PF_SimulateObserverParametric(StimLevels,params,...
%       OutOfNum,SDTfunc,M)
%
%returns something like:
%
%   NumPos =
% 
%     50    68    78    83    95   100
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

function NumPos = PAL_SDT_PF_SimulateObserverParametric(StimLevels,params,OutOfNum,SDTfunc,M)

nnums=length(StimLevels(1,:));
NumPos = zeros(1,nnums);

DP=(params(1).*StimLevels).^params(2);

if isempty(M)
    PC = SDTfunc(DP);
else
    PC = SDTfunc(DP,M);
end

for i=1:nnums
    Pos = rand(OutOfNum(i),1);
    Pos(Pos < PC(i)) =1;
    Pos(Pos ~= 1) = 0;
    NumPos(i) = sum(Pos);
end