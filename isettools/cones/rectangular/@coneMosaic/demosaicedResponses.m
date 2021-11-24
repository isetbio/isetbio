function varargout = demosaicedResponses(obj, varargin)
% Return demosaiced response maps (absorptions/currents) from a coneMosaic
%
% Syntax:
%   varargout = demoisaicedResponses(obj, [varargin])
%
% Description:
%    Return the demosaiced response maps (from absorptions or currents)
%    from a coneMosaic object.
%
% Inputs:
%    obj       - The cone mosaic object
%
% Outputs:
%    varargout - Array of varying size containing relevant output.
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% Notes:
%    * [Note: DHB - SOME HEADER COMMENTS HERE, PARTICULARLY ON THE INPUT
%      OPTIONS AND ALGORITHM, WOULD BE VERY HELPFUL.]
%

% History:
%    xx/xx/16  NPC  ISETBIO Team, 2016
%    02/23/18  jnm  Formatting

    if (isempty(varargin))
        responses = obj.absorptions;
    else
        responses = varargin{1};
    end

    coneMosaics = {'Lcone', 'Mcone', 'Scone'};
    coneIndices = containers.Map();
    coneIndices('null') = find(obj.pattern == 1);
    for coneMosaicIndex = 1:numel(coneMosaics)
        coneIndices(coneMosaics{coneMosaicIndex}) = ...
            find(obj.pattern == coneMosaicIndex + 1);
    end

    nFrames = size(responses, 3);
    demosaicedResponses = zeros([size(responses, 1), ...
        size(responses, 2), 3, nFrames]);

    % Set up interpolation mesh
    [XX, YY] = meshgrid(1:size(responses, 2), 1:size(responses, 1));

    % Interplate each frame
    for coneMosaicIndex = 1:numel(coneMosaics)
       if (numel(coneIndices(coneMosaics{coneMosaicIndex})) == 0)
           demosaicedResponses(:, :, coneMosaicIndex, 1:nFrames) = 0;
       else
         [coneRows, coneCols] = ind2sub([obj.rows obj.cols], ...
             coneIndices(coneMosaics{coneMosaicIndex}));
         for frameIndex = 1:nFrames
             responsesFrame = squeeze(responses(:, :, frameIndex));
             F = scatteredInterpolant(coneCols, coneRows, ...
                 responsesFrame(coneIndices(...
                 coneMosaics{coneMosaicIndex})), 'linear');
             demosaicedResponses(:, :, coneMosaicIndex, frameIndex) = ...
                 F(XX, YY);
             % demosaicedResponses(:, :, coneMosaicIndex, frameIndex) = ...
             %      griddata(coneCols, coneRows, ...
             %      responsesFrame(coneIndices(...
             %      coneMosaics{coneMosaicIndex})), XX, YY, 'linear');
         end % frameIndex
       end
    end % coneMosaicIndex

    % Return demosaiced maps
    varargout{1} = demosaicedResponses;

    % Return sRGB if we are operating on absorptions
    if (isempty(varargin))   
        sRGB = demosaicedResponses * 0;

        % Compute quantal efficiencies at the cornea
        lens = Lens;
        lensTransmittance = lens.transmittance;
        cornealQuantalEfficiencies = bsxfun(@times, obj.qe, ...
            lensTransmittance);

        % Find scaling so that conversion to sRGB is correct. That
        % conversion assumes that LMS spectral sensitivities are scaled to
        % peak of 1 in energy units. Here we have real quantal efficiences,
        % so that our LMS values are scaled differently relative to one
        % another. This will produce problems if we just use the standard
        % transformation matrix.
        %
        % Getting this totally right is a bit involved, particularly if we
        % want to get the XYZ values in standard units of Y being in cd/m2
        % or some such. Here we'll only fuss a little, and get the relative
        % scaling right.
        cornealQuantalEfficienciesEnergyUnits = EnergyToQuanta(...
            obj.wave(:), cornealQuantalEfficiencies);
        scale = max(cornealQuantalEfficienciesEnergyUnits);
        scale = scale / max(scale);
        scale = reshape(scale, [1 1 3]);

        % check number of non-zero entries to determine the mosaic is
        % trichromatic, or abnormal
        colorBlindType = find(~obj.spatialDensity(2:4));

        for frameIndex = 1:nFrames
            % scale LMS to match stockman normalized spectral qe. We need
            % to do this because all xyz2lms and colorblind transform
            % routine are based on stockman normalized spectral qe
            LMS = bsxfun(@rdivide, ...
                demosaicedResponses(:, :, :, frameIndex), scale);
            LMS = lms2lmsDichromat(LMS, colorBlindType, 'linear');
            sRGB(:, :, :, frameIndex) = lms2srgb(LMS);
        end

        % Return sRGB rendition
        varargout{2} = sRGB;
    end
end
