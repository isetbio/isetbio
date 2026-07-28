function valueOut = eccentricityInMacaqueRetina(conversionDirection, valueIn)

	switch (conversionDirection)
		case 'MMsToDegs'
			% valueIn is specified in MMs -> valueOut will be in degs
			valueOut = valueIn * 1e3 / RGCmodels.CronerKaplan.constants.micronsPerDegreeInMacaqueRetina;

		case 'MicronsToDegs'
			% valueIn is specified in microns -> valueOut will be in degs
			valueOut = valueIn / RGCmodels.CronerKaplan.constants.micronsPerDegreeInMacaqueRetina;

		case 'DegsToMMs'
			% valueIn in specified in visual degs -> valueOut will be in MMs
			valueOut = 1e-3 * valueIn * RGCmodels.CronerKaplan.constants.micronsPerDegreeInMacaqueRetina;

		case 'DegsToMicrons'
			% valueIn in specified in visual degs -> valueOut will be in MMs
			valueOut = valueIn * RGCmodels.CronerKaplan.constants.micronsPerDegreeInMacaqueRetina;

		otherwise
			fprintf(2, 'Valid conversions are {''MMsToDegs'', ''MicronsToDegs'', ''DegsToMMs'', and ''DegsToMicrons''}.')
			error('Unknown conversion direction: ''%s''.', conversionDirection);
	end
end