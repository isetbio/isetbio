
function fracVar = calculateFractionalVariance(innerRetinaPSTH, innerRetinaRecordedPSTH, stimulusTestI)
% Calculate the fractional variance for the predictions in Fig. 2.
% 
% Inputs: simulated and recorded PSTHs, stimulusTestI
% 
% 5/2016 JRG (c) isetbio team

nCells = size(innerRetinaPSTH,2);
fracVar = zeros(nCells,1);

minlen = min([length(innerRetinaPSTH{1}) length(innerRetinaRecordedPSTH{1}) ]);

for i = 1:nCells

switch stimulusTestI
    case 1
        rsim = innerRetinaPSTH{i}(600+(1:minlen-1200));
    case 2
        rsim = innerRetinaPSTH{i}(1200+(1:minlen-1200));
end

rrec = innerRetinaRecordedPSTH{i}(1:minlen-1200);
fracVar(i) = 1 - sum((rsim-rrec).^2)/sum((rrec-mean(rrec)).^2);


end

% % Calculate inter-trial reliability J
% 
% innerRetinaRecordedOdd = irPhys(os1, params);
% innerRetinaRecordedOdd = irSet(innerRetinaRecordedOdd,'numberTrials',nTrials);
% 
% iTctr = 0;
% for cellind = 1:length(xval_mosaic)
%     for iTrial = 1:2:nTrials
%         iTctr = iTctr+1;
%         recorded_spiketimes{1,cellind,iTctr} = find(xval_mosaic{cellind}.rasters.recorded(iTrial,:)>0);
%     end
% end
% 
% innerRetinaRecordedOdd = irSet(innerRetinaRecordedOdd,'recordedSpikes',recorded_spiketimes);
% innerRetinaRecordedOddPSTH = mosaicGet(innerRetinaRecordedOdd.mosaic{1},'responsePsth');
% % % % % % 
% innerRetinaRecordedEven = irPhys(os1, params);
% innerRetinaRecordedEven = irSet(innerRetinaRecordedEven,'numberTrials',nTrials);
% 
% iTctr = 0;
% for cellind = 1:length(xval_mosaic)
%     for iTrial = 2:2:nTrials
%         iTctr = iTctr+1;
%         recorded_spiketimes{1,cellind,iTctr} = find(xval_mosaic{cellind}.rasters.recorded(iTrial,:)>0);
%     end
% end
% 
% innerRetinaRecordedEven = irSet(innerRetinaRecordedEven,'recordedSpikes',recorded_spiketimes);
% innerRetinaRecordedEvenPSTH = mosaicGet(innerRetinaRecordedEven.mosaic{1},'responsePsth');
% for i = 1:68
%     J(i) = 1 - sum((innerRetinaRecordedEvenPSTH{i}-innerRetinaRecordedOddPSTH{i}).^2)/sum((innerRetinaRecordedEvenPSTH{i}-mean(innerRetinaRecordedEvenPSTH{i})).^2)
% end