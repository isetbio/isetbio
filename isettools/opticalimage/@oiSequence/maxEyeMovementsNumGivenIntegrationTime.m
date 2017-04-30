function maxEyeMovementsNum = maxEyeMovementsNumGivenIntegrationTime(obj,integrationTime)
    % Generate eye movement sequence for all oi's
    if (numel(obj.timeAxis) == 1)
        stimulusSamplingInterval  = integrationTime;
    else
        stimulusSamplingInterval = obj.timeAxis(2)-obj.timeAxis(1);
    end
      
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/integrationTime;
    maxEyeMovementsNum = round(eyeMovementsNumPerOpticalImage*obj.length);

    if (maxEyeMovementsNum < 1)
        error('Less than 1 eye movement per oiSequence !!! \nStimulus sampling interval:%g ms Cone mosaic integration time: %g ms\n', 1000*stimulusSamplingInterval, 1000*integrationTime);
    else 
        fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', maxEyeMovementsNum, eyeMovementsNumPerOpticalImage);
    end 
end