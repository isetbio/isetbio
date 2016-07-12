%
%PAL_PFML_GroupTrialsbyX    Combines like trials (i.e., those with 
%   identical stimulus levels) in psychometric function fitting setting.
%
%syntax: [StimLevelsG NumPosG OutOfNumG] = ...
%   PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum)
%
%example:
%
%   StimLevels = [2 4 3 1 3];   %stimulus levels (3 is listed twice)
%   NumPos =     [1 2 3 4 2];   %number of trials with positive response
%                               %at each level (total of 5 at StimLevels ==
%                               %3)
%   OutOfNum =   [5 10 10 5 5]; %number of trials at each level (total of
%                               %15 at StimLevels == 3)
%
%   [StimLevelsG NumPosG OutOfNumG] = ...
%      PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum) returns:
%
%   StimLevelsG =  1     2     3     4
%
%   NumPosG =      4     1     5     2
%
%   OutOfNumG =    5     5    15    10
%
%Input variables may also be matrices (for example when results are those 
%   of a multi-condition experiment). In such cases trials will be combined
%   within each row of the matrices separately. Example:
%
%   StimLevels =    [0  1  2  3  1;
%                    7  6  5  4  7];
%
%   NumPos =        [0  5  5  5  5;
%                    10 10 10 10 10];
%
%   OutOfNum =      [0  10 10 10 10;
%                    20 20 20 20 20];
%
%   [StimLevelsG NumPosG OutOfNumG] = ...
%       PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum) returns:
%
%   StimLevelsG =
%    1     2     3     0
%    4     5     6     7
%
%   NumPosG =
%   10     5     5     0
%   10    10    10    20
%
%   OutOfNumG =
%   20    10    10     0
%   20    20    20    40
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.2.0, 1.6.3 (see History.m)

function [StimLevelsG, NumPosG, OutOfNumG] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum)

for row = 1:size(StimLevels,1)

    StimLevelsRow = squeeze(StimLevels(row,:));
    NumPosRow = squeeze(NumPos(row,:));
    OutOfNumRow = squeeze(OutOfNum(row,:));
    
    StimLevelsRow = StimLevelsRow(OutOfNumRow ~= 0);
    NumPosRow = NumPosRow(OutOfNumRow ~= 0);
    OutOfNumRow = OutOfNumRow(OutOfNumRow ~= 0);
    NumPosGRow = [];
    OutOfNumGRow = [];
    xGRow = [];
    
    NumLevels = 1;

    while length(StimLevelsRow) ~= 0

        x = StimLevelsRow(1)*ones(size(StimLevelsRow));
        xGRow(NumLevels) = StimLevelsRow(1);
        match = (StimLevelsRow == x);
        OutOfNumGRow(NumLevels) = sum(OutOfNumRow(match));
        NumPosGRow(NumLevels) = sum(NumPosRow(match));
        OutOfNumRow = OutOfNumRow(match == 0);
        NumPosRow = NumPosRow(match == 0);
        StimLevelsRow = StimLevelsRow(match == 0);
        NumLevels = NumLevels+1;

    end

    [xG I] = sortrows(xGRow');
    NumPosG(row,1:length(NumPosGRow)) = NumPosGRow(I);
    OutOfNumG(row,1:length(OutOfNumGRow)) = OutOfNumGRow(I);
    StimLevelsG(row,1:length(xG)) = xG';

end