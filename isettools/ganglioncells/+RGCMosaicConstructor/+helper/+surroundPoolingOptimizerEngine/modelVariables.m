function modelVariablesStruct = modelVariables(modelConstants, ...
  anatomicalRcDegs, poolingModel, userSuppliedInitialValues)

  names = {}; plotScaling = {}; 
  initialValues = []; lowerBounds = []; upperBounds = [];

  names{numel(names)+1} = 'Kc';
  plotScaling{numel(plotScaling)+1} = 'log';
  initialValues(numel(initialValues)+1) = 1;
  lowerBounds(numel(lowerBounds)+1) = 1;
  upperBounds(numel(upperBounds)+1) = 1;

  names{numel(names)+1} = 'RwDegs';
  RnarrowToRwideRatiosExampleH1Cells = RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.RnarrowToRwideRatios;
  centerConesNum = numel(find(modelConstants.cachedData.centerConeWeights >= mRGCMosaic.sensitivityAtPointOfOverlap));
  RwDegsInitial = 2.00 * 1/mean(RnarrowToRwideRatiosExampleH1Cells) * anatomicalRcDegs * sqrt(2.3) * sqrt(centerConesNum);
  plotScaling{numel(plotScaling)+1} = 'linear';

  initialValues(numel(initialValues)+1) = RwDegsInitial;
  lowerBounds(numel(lowerBounds)+1) = max([modelConstants.cachedData.minSurroundRadiusDegs  0.0001*RwDegsInitial]);
  upperBounds(numel(upperBounds)+1) = max([modelConstants.cachedData.maxSurroundSupportDegs RwDegsInitial]);

  names{numel(names)+1} = 'KsKcRatio';
  plotScaling{numel(plotScaling)+1} = 'log';
  initialValues(numel(initialValues)+1) = 0.06;
  lowerBounds(numel(lowerBounds)+1) = 0.005;
  upperBounds(numel(upperBounds)+1) = 1.0;

  switch (poolingModel.name)
    case 'PackerDacey2002H1free'
    
      names{numel(names)+1} = 'VnVwRatio';
      plotScaling{numel(plotScaling)+1} = 'log';
      initialValues(numel(initialValues)+1) = mean(poolingModel.H1parameterTolerances.VnarrowToVwideRatio);
      if (numel(poolingModel.H1parameterTolerances.VnarrowToVwideRatio) == 2)
        lowerBounds(numel(lowerBounds)+1) = poolingModel.H1parameterTolerances.VnarrowToVwideRatio(1);
        upperBounds(numel(upperBounds)+1) = poolingModel.H1parameterTolerances.VnarrowToVwideRatio(2);
      else
        error('poolingModel.params.H1parameterTolerances.VnarrowToVwideRatio must be either a 2-element vector ');
      end

      names{numel(names)+1} = 'RnRwRatio';
      plotScaling{numel(plotScaling)+1} = 'log';
      initialValues(numel(initialValues)+1) = mean(poolingModel.H1parameterTolerances.RnarrowToRwideRatio);
      if (numel(poolingModel.H1parameterTolerances.RnarrowToRwideRatio) == 2)
        lowerBounds(numel(lowerBounds)+1) = poolingModel.H1parameterTolerances.RnarrowToRwideRatio(1);
        upperBounds(numel(upperBounds)+1) = poolingModel.H1parameterTolerances.RnarrowToRwideRatio(2);
      else
        error('poolingModel.params.H1parameterTolerances.RnarrowToRwideRatio must be a  2-element vector');
      end

      if (~isempty(userSuppliedInitialValues))
        initialValues = updateInitialValuesUsingUserSuppliedInitialValues(...
          poolingModel, names, initialValues, userSuppliedInitialValues);
      else
        fprintf('No user supplied initial values for model variables. Using defaults.\n');
      end

    case 'PackerDacey2002H1FixedCellIndex'

      fixedH1index = poolingModel.fixedH1CellIndex;

      names{numel(names)+1} = 'VnVwRatio';
      plotScaling{numel(plotScaling)+1} = 'log';
      initialValues(numel(initialValues)+1) = RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.VnarrowToVwideRatios(fixedH1index);
      if (numel(poolingModel.H1parameterTolerances.VnarrowToVwideRatio) == 1)
        lowerBounds(numel(lowerBounds)+1) = max([0.01 initialValues(end)-poolingModel.H1parameterTolerances.VnarrowToVwideRatio]);
        upperBounds(numel(upperBounds)+1) = initialValues(end)+poolingModel.H1parameterTolerances.VnarrowToVwideRatio;
      elseif (numel(poolingModel.params.H1parameterTolerances.VnarrowToVwideRatio) == 2)
        lowerBounds(numel(lowerBounds)+1) = max([0.01 initialValues(end)-poolingModel.H1parameterTolerances.VnarrowToVwideRatio(1)]);
        upperBounds(numel(upperBounds)+1) = initialValues(end)+poolingModel.H1parameterTolerances.VnarrowToVwideRatio(2);
      else
        error('poolingModel.H1parameterTolerances.VnarrowToVwideRatio must be either a 1- or a 2-element vector');
      end


      names{numel(names)+1} = 'RnRwRatio';
      plotScaling{numel(plotScaling)+1} = 'log';
      initialValues(numel(initialValues)+1) = RGCMosaicConstructor.constants.PackerDaceyParamsForH1CellIndices.RnarrowToRwideRatios(fixedH1index);
      if (numel(poolingModel.H1parameterTolerances.RnarrowToRwideRatio) == 1)
        lowerBounds(numel(lowerBounds)+1) = max([0.01 initialValues(end)-poolingModel.H1parameterTolerances.RnarrowToRwideRatio]);
        upperBounds(numel(upperBounds)+1) = initialValues(end)+poolingModel.H1parameterTolerances.RnarrowToRwideRatio;
      elseif (numel(poolingModel.params.H1parameterTolerances.RnarrowToRwideRatio) == 2)
        lowerBounds(numel(lowerBounds)+1) = max([0.01 initialValues(end)-poolingModel.H1parameterTolerances.RnarrowToRwideRatio(1)]);
        upperBounds(numel(upperBounds)+1) = initialValues(end)+poolingModel.H1parameterTolerances.RnarrowToRwideRatio(2);
      else
        error('poolingModel.H1parameterTolerances.RnarrowToRwideRatio must be either a 1- or a 2-element vector');
      end

      if (~isempty(userSuppliedInitialValues))
        initialValues = updateInitialValuesUsingUserSuppliedInitialValues(...
          poolingModel, names, initialValues, userSuppliedInitialValues);
      else
        fprintf('No user supplied initial values for model variables. Using defaults.\n');
      end

    otherwise
      error('Model ''%s'' not implemented.', poolingModel.name);
  end % switch

  modelVariablesStruct.names = names;
  modelVariablesStruct.scaling = plotScaling;
  modelVariablesStruct.initialValues = initialValues;
  modelVariablesStruct.lowerBounds = lowerBounds;
  modelVariablesStruct.upperBounds = upperBounds;
end 

function initialValues = updateInitialValuesUsingUserSuppliedInitialValues(...
          poolingModel, names, initialValues, userSuppliedInitialValues)

  % Verify that the userSuppliedInitialValuesForModelVariables are for this model
  assert(isstruct(userSuppliedInitialValues), ...
    'userSuppliedInitialValuesForModelVariables is not a struct.');
  assert(isfield(userSuppliedInitialValues, 'poolingModel'), ...
    'userSuppliedInitialValuesForModelVariables does not contain a ''poolingModel'' field.');
  assert(strcmp(userSuppliedInitialValues.poolingModel.name, poolingModel.name), ...
    'userSuppliedInitialValuesForModelVariables is for a different model (''%s'') than the current one (''%s'').', userSuppliedInitialValues.poolingModel.name, poolingModel.name);
  assert(isequal(userSuppliedInitialValues.poolingModel.weightsComputeHandle, poolingModel.weightsComputeHandle), ...
    'the @weightsComputeFunctionHandle in the  userSuppliedInitialValuesForModelVariables is different than the function handle for the current model');
  assert(isfield(userSuppliedInitialValues, 'optimizedValues'), ...
    'userSuppliedInitialValuesForModelVariables does not contain an ''optimizedValues'' field.');

  % All OK. Overwrite the current models variables with the supplied values
  variablesNames = fieldnames(userSuppliedInitialValues.optimizedValues);
  for iVar = 1:numel(variablesNames)
    vName = variablesNames{iVar};
    idx = find(strcmp(names, vName));
    if (isempty(idx))
      error('User supplied variable name ''%s'' is not part of the current model''s variables.', vName);
    end
    newInitialValue = userSuppliedInitialValues.optimizedValues.(vName);
    fprintf('Replacing initial value for ''%s'' (%f) with %f (user-supplied)\n', vName, initialValues(idx), newInitialValue);
    initialValues(idx) = newInitialValue;
  end % iVar
end
