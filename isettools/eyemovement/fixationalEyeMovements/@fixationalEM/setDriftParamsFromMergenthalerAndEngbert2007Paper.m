function setDriftParamsFromMergenthalerAndEngbert2007Paper(obj)
    % Control gain (gamma in Mergenthaler&Engbert, 2007)
    obj.controlGamma = 0.25;
    
    % Control noise, ksi, with zero MEAN and sigma STD
    ksi = struct('mean', 0, 'sigma', 0.075);
    obj.controlNoiseMean = ksi.mean;
    obj.controlNoiseStd = ksi.sigma; 
    
    % Position noise, eta with zero MEAN and rho STD
    eta = struct('mean', 0, 'rho', 0.35);
    obj.positionNoiseMean = eta.mean;
    obj.positionNoiseStd = eta.rho;    
    
    %Feedback gain (lambda in Mergenthaler&Engbert, 2007)
    obj.feedbackGain = 0.15;             
     
    % Feedback steepness (epsilon in Mergenthaler&Engbert, 2007)
    obj.feedbackSteepness = 1.1;          
    
    % Feedback delay for x- and y-pos 
    obj.feedbackXposDelaySeconds = 70/1000; % tau_hor in Mergenthaler&Engbert, 2007
    obj.feedbackYposDelaySeconds = 40/1000; % tau_vert in Mergenthaler&Engbert, 2007
end