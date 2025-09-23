function [eccDegsOut, eccMMsOut] = eccentricityDegsBetweenSpecies(conversionDirection, eccDegsIn)

	switch (conversionDirection)
		case 'HumanRetinaToMacaqueRetina'
			% Convert eccDegs to eccMMs in human retina
			eccMMsIn = RGCmodels.Watson.convert.rhoDegsToMMs(eccDegsIn);

			% Assume retinal MMs are equivalent in the 2 species
			eccMMsOut = eccMMsIn;

			% Convert eccMMs to eccDegs in macaque retina
			eccDegsOut = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('MMsToDegs', eccMMsOut);

		case 'MacaqueRetinaToHumanRetina'
			% Convert eccDegs to eccMMs in macaque retina
			eccMMsIn = RGCMosaicConstructor.helper.convert.eccentricityInMacaqueRetina('DegsToMMs', eccDegsIn);

			% Assume retinal MMs are equivalent in the 2 species
			eccMMsOut = eccMMsIn;

			% Convert eccMMs to eccDegs in human retina
			eccDegsOut = RGCmodels.Watson.convert.rhoMMsToDegs(eccMMsOut);

		otherwise
			fprintf(2, 'Valid conversions are {''HumanRetinaToMacaqueRetina'' and ''MacaqueRetinaToHumanRetina''}.');
			error('Unknown conversion direction: ''%s''.', direction);
	end % switch
end
