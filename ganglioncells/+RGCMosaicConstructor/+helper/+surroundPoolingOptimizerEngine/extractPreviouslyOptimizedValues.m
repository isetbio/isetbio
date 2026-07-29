function optimizedValuesStruct = extractPreviouslyOptimizedValues(optimizationResults)
	% Assemble the optimizedValuesStruct
	optimizedValuesStruct = struct();

	for iVar = 1:numel(optimizationResults.modelVariables.names)
	    theVariableName = optimizationResults.modelVariables.names{iVar};
	    optimizedValuesStruct.(theVariableName) = optimizationResults.modelVariables.finalValues(iVar);
	end % iVar
end
