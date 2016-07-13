%
%PAL_SDT_Summ_MultiplePFML_negLL     (negative) Log Likelihood associated with fit 
%   of summation psychometric function 
%
%Internal Function
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

function negLL = PAL_SDT_Summ_MultiplePFML_negLL(params,StimLevels,NumPos,OutOfNum,SummFunc,M,Q)

[nrows, ncols, nnums]=size(StimLevels);
PC=zeros(nrows,nnums);
StimVec=zeros(1,ncols); % vector for corresponding stimulus levels for summation equation

num=length(params);
gParams=params(1:num/2);
pParams=params(num/2+1:num);

for j=1:nrows
    for i=1:nnums
        for k=1:ncols
            StimVec(1,k)=StimLevels(j,k,i);
        end
        PC(j,i) = SummFunc(StimVec,gParams,pParams,M,Q);
    end
end


negLL = -sum(PAL_nansum(NumPos.*log(PC))+PAL_nansum((OutOfNum-NumPos).*log(1 - PC)));