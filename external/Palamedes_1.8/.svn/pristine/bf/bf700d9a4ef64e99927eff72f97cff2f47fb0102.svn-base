%
%PAL_MLDS_GenerateStimList   Generate stimulus sets of pairs, triads or 
%   quadruples for use in scaling experiment.
%
%   syntax: StimList = PAL_MLDS_GenerateStimList(ptq, NumLevels, MaxDiff, 
%        NumRep)
%
%   'ptq' should be either 2, 3, or 4 depending on whether one wishes to
%       generate pairs, triads or quadruples respectively.
%
%   'NumLevels' is the number of stimulus levels one wishes to utilize in
%       experiment.
%
%   'MaxDiff' precludes stimulus combinations that are 'too far apart'. 
%       Specifically, if 'ptq' is set to 2, all stimulus pairs i-j (i < j) 
%       for which (j - i) has a value not exceeding 'MaxDiff' will be 
%       included in stimulus set. If 'ptq' is set to 3, all triads i-j-k 
%       (i < j < k) for which |(k - j) - (j - i)| has a value not exceeding
%       'MaxDiff' will be included in stimulus set. If 'ptq' is set to 4, 
%       all quadruples i-j-k-l (i < j < k < l) for which 
%       |(l - k) - (j - i)| has a value not exceeding 'MaxDiff' will be 
%       included in stimulus set.
%
%    'NumRep' sets the number of times each possible pair/triad/quadruple 
%       will be included in the stimulus set.
%
%   Example:
%
%       StimList = PAL_MLDS_GenerateStimList(2, 4, 1, 2)
%
%       returns:
%
%       StimList =
%
%        1     2
%        2     3
%        3     4
%        1     2
%        2     3
%        3     4
%
%   PAL_randomizeArray may be used to obtain a randomized stimulus list,
%   e.g.:
%
%   StimList = PAL_randomizeArray(PAL_MLDS_GenerateStimList(2, 4, 1, 2))
%
%   will return the same list as above but now in a random order.
%
%Introduced: Palamedes version 1.0.0 (NP)

function StimList = PAL_MLDS_GenerateStimList(ptq, NumLevels, MaxDiff, NumRep)

count = 0;

if ptq == 2
    for i = 1:NumLevels-1
        for j = i+1:min(i+MaxDiff,NumLevels)
            count = count+1;
            StimList(count,:) = [i j];
        end
    end
end
if ptq == 3
    for i = 1:NumLevels-2
        for j = i+1:NumLevels-1
            for k = j+1:NumLevels
                count = count+1;
                StimList(count,:) = [i j k];
            end
        end
    end
    StimList = StimList(abs(StimList(:,3)-2*StimList(:,2)+StimList(:,1))<=MaxDiff,:);
end
if ptq == 4
    for i = 1:NumLevels-3
        for j = i+1:NumLevels-2
            for k = j+1:NumLevels-1
                for l = k+1:NumLevels
                    count = count+1;
                    StimList(count,:) = [i j k l];
                end
            end
        end
    end
    StimList = StimList(abs(StimList(:,4)-StimList(:,3)-StimList(:,2)+StimList(:,1))<=MaxDiff,:);
end
StimList = repmat(StimList,[NumRep 1]);