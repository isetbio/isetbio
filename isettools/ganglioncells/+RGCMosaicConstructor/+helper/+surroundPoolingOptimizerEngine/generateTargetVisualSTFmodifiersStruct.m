function targetVisualSTFmodifierStruct = generateTargetVisualSTFmodifiersStruct(targetVisualSTFdescriptor)

switch (targetVisualSTFdescriptor)
	case 'default'
	% Default arget visual STF modifiers
	targetVisualSTFmodifierStruct = struct(...
        'surroundToCenterRcRatioMultiplier', 1.0, ...
        'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.0);

    case 'x1.5 RsRcRatio'
        % Target visual STF modifiers: 1.5x multiplier for surround to center characteristic radius ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.5, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.0);

    case 'x1.3 RsRcRatio'
    	% Target visual STF modifiers: 1.3 multiplier for surround to center characteristic radius ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.3, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.0);

	case 'x0.75 RsRcRatio'
    	% Target visual STF modifiers: 1.3 multiplier for surround to center characteristic radius ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 0.75, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.0);

	case 'higher intSCratio'
		% Target visual STF modifiers: 1.3 multiplier for surround to center integrated sensitivity ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.0, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.3);

	case 'very high intSCratio'
		% Target visual STF modifiers: 1.7 multiplier for surround to center integrated sensitivity ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.0, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 1.7);

	case 'ultra high intSCratio'
		% Target visual STF modifiers: 1.7 multiplier for surround to center integrated sensitivity ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.0, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 2.0);

	case 'lower intSCratio'
		% Target visual STF modifiers: 0.75 multiplier for surround to center integrated sensitivity ratio
		targetVisualSTFmodifierStruct = struct(...
        	'surroundToCenterRcRatioMultiplier', 1.0, ...
        	'surroundToCenterIntegratedSensitivityRatioMultiplier', 0.75);

    otherwise
        error('Unknown targetVisualSTF: ''%s'' !!', targetVisualSTF);
end % switch (targetVisualSTF)

end