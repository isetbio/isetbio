function unpausingExecution(theMessage)
%
% RGCMosaicConstructor.helper.queryUserFor.unpausingExecution(theMessage)
%
	% Message user
	if (nargin == 0)
		fprintf('\nPaused. Hit enter to continue...\n');
	else
		fprintf('\n%s. Hit enter to continue...\n', theMessage);
	end

	% And pause
	pause;

end