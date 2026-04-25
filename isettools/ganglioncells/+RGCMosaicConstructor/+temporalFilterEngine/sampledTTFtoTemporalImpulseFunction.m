%
% RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
%    theTTF, theFrequencySupport)
%

function theImpulseResponseFunctionStruct = sampledTTFtoTemporalImpulseFunction(theTTF, theFrequencySupportHz, varargin)

% Compute the impulse response from a sampled transfer function (one-sided or two-sided, real or complex).
%
% RESOLUTION OPTIONS
%   Zero-padding the TIME domain  ('zeropad')   -> finer frequency axis
%   Zero-padding the FREQ domain  ('upsample')  -> finer time axis
%
%   These are the exact duals of each other:
%     zeropad  N->N*zp  : dt unchanged,  df -> df/zp   (more freq points)
%     upsample N->N*us  : df unchanged,  dt -> dt/us   (more time points)
%
%   upsample inserts zeros between existing spectral bins (equivalent to
%   sinc-interpolating the impulse response in time).
%
%
% INPUTS
%   theTTF   - Sampled transfer function, length M (complex or magnitude-only)
%   theFrequencySupportHz   - Frequency axis (Hz), same length as H
%         Accepted formats:
%           one-sided : f = [0, df, ..., fs/2]             length N/2+1
%           two-sided : f = [0, df, ..., -df]              length N
%           centered  : f = [-fs/2, ..., 0, ..., fs/2]    length N
%
% NAME-VALUE OPTIONS
%   'format'     - 'auto'|'onesided'|'twosided'|'centered'  (default:'auto')
%   'window'     - 'none'|'hann'|'hamming'|'blackman'|'kaiser'|'tukey'
%                                                            (default:'none')
%   'beta'       - Kaiser beta                               (default: 8)
%   'alpha'      - Tukey taper fraction                      (default: 0.5)
%   'causal'     - true : output starts at t=0              (default: false)
%   'normalize'  - Normalize peak of |h| to 1               (default: false)
%   'trim'       - Trim h to where envelope > trim_thr       (default: false)
%   'trim_thr'   - Amplitude threshold for trimming          (default: 0.01)
%   'fs'         - Sample rate (Hz). Inferred from f if absent.
%
%   ── Frequency-domain zero-padding (finer TIME axis) 
%   'upsample'   - Integer factor >= 1.
%                  Inserts zeros into the SPECTRUM, extending the effective
%                  sample rate to fs*upsample.
%                  dt_new = 1/(fs*upsample)
%                  The output h and t reflect the upsampled resolution.
%                                                            (default: 1)
%   'upsample_method'
%                - 'zeroinsert' : insert spectral zeros at high freqs
%                                 (band-limited / sinc interpolation)
%                  'interp'     : MATLAB interp() on the impulse response
%                                 (faster for large factors, less precise)
%                                                         (default:'zeroinsert')
%
%   ── Time-domain zero-padding (finer FREQUENCY axis) 
%   'zeropad'    - Integer factor >= 1.
%                  Zero-pads the impulse response before re-FFT.
%                  
%

    % Parse input
    p = inputParser;
    p.addRequired('theTTF');
    p.addRequired('theFrequencySupportHz');
    p.addParameter('format',          'auto');
    p.addParameter('window',          'none');
    p.addParameter('beta',            8);
    p.addParameter('alpha',           0.5);
    p.addParameter('causal',          true);
    p.addParameter('normalize',       false);
    p.addParameter('trim',            false);
    p.addParameter('trim_thr',        0.01);
    p.addParameter('samplingFrequency',              []);
    p.addParameter('upsample',        4);
    p.addParameter('upsample_method', 'zeroinsert');
    p.addParameter('zeropad',         1);
    p.addParameter('plot',            false);
    p.addParameter('beVerbose', false);
    p.parse(theTTF, theFrequencySupportHz, varargin{:});
    options = p.Results;
    
    assert(options.upsample >= 1 && floor(options.upsample)==options.upsample, ...
        '''upsample'' must be a positive integer.');
    assert(options.zeropad  >= 1 && floor(options.zeropad) ==options.zeropad,  ...
        '''zeropad'' must be a positive integer.');

    theTTF = theTTF(:);
    theFrequencySupportHz = theFrequencySupportHz(:);
    M = length(theTTF);
    assert(length(theFrequencySupportHz)==M,'theTTF and theFrequencySupportHz must have the same length.');

    % Data format
    fmt = options.format;
    if strcmp(fmt,'auto')
        if theFrequencySupportHz(1) >= 0 && theFrequencySupportHz(end) > 0 && all(diff(theFrequencySupportHz) > 0)
            fmt = 'onesided';
        elseif theFrequencySupportHz(1) < 0 || (theFrequencySupportHz(1)==0 && any(theFrequencySupportHz < 0))
            fmt = 'centered';
        else
            fmt = 'twosided';
        end
    end

    % Sample rate
    switch fmt
        case 'onesided'
            df = theFrequencySupportHz(2) - theFrequencySupportHz(1);
            samplingFrequency = 2 * theFrequencySupportHz(end);
            N  = 2 * (M - 1);
        case {'centered','twosided'}
            df = theFrequencySupportHz(2) - theFrequencySupportHz(1);
            N  = M;
            samplingFrequency = N * df;
    end
    if (~isempty(options.samplingFrequency))
        samplingFrequency = options.samplingFrequency;
        deltaFrequency = samplingFrequency/N;
    end
    sampleIntervalSeconds = 1/samplingFrequency;

    if (options.beVerbose)
        fprintf('──────────────────────────────────────────────\n');
        fprintf('Format    : %s\n', fmt);
        fprintf('N         : %d   delta-frequency = %.5fHz   samplingFrequency = %.2fHz   dt = %.5f seconds\n', N, deltaFrequency, samplingFrequency, sampleIntervalSeconds );
        if options.upsample > 1
            fprintf('Upsample  : %dx  [%s]  ->  sampling interval = %.6f seconds   upsampled frequency = %.2f Hz\n', ...
                options.upsample, options.upsample_method, sampleIntervalSeconds/options.upsample, samplingFrequency * options.upsample);
        end
        if options.zeropad > 1
            fprintf('Zero-pad  : %d x  ->  df_hires = %.5f Hz\n', ...
                options.zeropad, deltaFrequency/options.zeropad);
        end
        fprintf('──────────────────────────────────────────────\n');
    end

    % Spectral window 
    win = generateWindow(N, options.window, options.beta, options.alpha);

    % Build two-sided FFT spectrum
    switch fmt
        case 'onesided'
            theTTF_half          = theTTF;
            theTTF_half(2:end-1) = theTTF(2:end-1) / 2;
            mirror               = conj(theTTF_half(end-1:-1:2));
            theTTF_std           = [theTTF_half; mirror];
            f_std                = [theFrequencySupportHz; -(theFrequencySupportHz(end-1:-1:2))];
            win                  = ifftshift(win);
        case 'centered'
            theTTF_std = ifftshift(theTTF);
            win        = ifftshift(win);
            f_std      = ifftshift(theFrequencySupportHz);
        case 'twosided'
            theTTF_std = theTTF;
            f_std      = theFrequencySupportHz;
    end

    % Make the FFT spectrum symmetric
    theTTF_herm           = theTTF_std;
    theTTF_herm(2:N/2)    = theTTF_std(2:N/2);
    theTTF_herm(N/2+2:N)  = conj(theTTF_std(N/2:-1:2));
    theTTF_herm(1)        = real(theTTF_std(1));
    theTTF_herm(N/2+1)    = real(theTTF_std(N/2+1));
    
    
    windowedSpectrum = theTTF_herm .* win;

    if (options.plot)
        figure(33)
        subplot(2,1,1);
        plot(f_std, abs(theTTF_std), 'ko');

        subplot(2,1,2)
        plot(abs(theTTF_herm), 'ko');
        hold on;
        plot(win, 'k-')
        pause
    end

    % UPSAMPLE: zero-insert in spectrum -> finer time axis ─────────────────────
    %
    % Upsampling by factor U means the new sample rate is fs*U and the new
    % FFT length is N*U.  The spectrum of the upsampled signal is obtained by
    % inserting (U-1)*N/2 zeros between the positive and negative frequency
    % halves (i.e., extending the band from fs/2 to fs*U/2 with zeros).
    %
    % Standard FFT layout:  [DC | pos 1..N/2-1 | Nyq | neg N/2+1..N-1]
    %                         ^                   ^
    %                         1                  N/2+1
    %
    % After inserting zeros at high frequencies:
    %   [DC | pos 1..N/2-1 | zeros (U-1)*N | Nyq | neg N/2+1..N-1]
    % New length: N*U  (the Nyquist bin is kept at the center)

    options.upsample = options.upsample;
    N_up = N * options.upsample;
    fs_up = samplingFrequency * options.upsample;

    if (options.upsample > 1) && (strcmpi(options.upsample_method,'zeroinsert'))
        % Number of zeros to insert between pos and neg halves
        n_zeros  = (options.upsample - 1) * N;
        upsampledSpectrum  = [windowedSpectrum(1:N/2);           ... DC through positive freqs
                              zeros(n_zeros, 1);       ... inserted high-freq zeros
                              windowedSpectrum(N/2+1:end)];       ... Nyquist through negative freqs
        % Scale to preserve amplitude (upsampling by options.upsample scales IFFT by options.upsample)
        upsampledSpectrum  = upsampledSpectrum * options.upsample;
    
    elseif (options.upsample > 1) && (strcmpi(options.upsample_method,'interp'))
        % Fast path: IFFT at native rate, then MATLAB interp()
        h_native  = real(ifft(windowedSpectrum));
        h_interp  = interp(h_native, options.upsample);     % lowpass interpolation
        h_interp  = h_interp(1:N_up);        % trim to exact length
        upsampledSpectrum = fft(h_interp);
    else
        upsampledSpectrum      = windowedSpectrum;                   % no upsampling
        N_up  = N;
        fs_up = fs;
    end

    % IFFT -> impulse response at upsampled rate 
    theImpulseResponseFunction = real(ifft(upsampledSpectrum));
    
    if (~options.causal)
        [~, idx] = max(abs(theImpulseResponseFunction));
        theImpulseResponseFunction = circshift(theImpulseResponseFunction, N_up/2 - idx + 1);
    end

    if (options.normalize)
        theImpulseResponseFunction = theImpulseResponseFunction / max(abs(theImpulseResponseFunction));
    end
    
    if (options.causal)
        temporalSupportSeconds = (0:N_up-1)' / fs_up;
    else
        temporalSupportSeconds = (-N_up/2 : N_up/2-1)' / fs_up;
    end

    % ZEROPAD: time-domain zero-padding -> finer frequency axis
    % This operates on the NATIVE (pre-upsample) impulse response for
    % spectral inspection; does not affect output.
    
    h_native_zp = real(ifft(windowedSpectrum));        % native rate, for zeropad branch
    if (~options.causal)
        [~, idx]    = max(abs(h_native_zp));
        h_native_zp = circshift(h_native_zp, N/2 - idx + 1);
    end

    N_pad   = N * options.zeropad;
    df_hi   = samplingFrequency / N_pad;

    if (options.zeropad) > 1
        if (options.causal)
            n_zeros  = N * (options.zeropad - 1);
            h_pad    = [h_native_zp(1:N/2); zeros(n_zeros,1); h_native_zp(N/2+1:end)];
        else
            h_pad    = [h_native_zp; zeros(N*(options.zeropad-1), 1)];
        end
        theTTF_hires  = fft(h_pad);
        f_hires  = (0 : N_pad/2)' * df_hi;
        theTTF_hires  = theTTF_hires(1:N_pad/2+1);
    else
        theTTF_hires  = theTTF_herm(1:N/2+1);
        f_hires  = (0:N/2)' * df;
        df_hi    = df;
    end

    % Trim
    theImpulseResponseEnvelope = abs(hilbert(theImpulseResponseFunction));
    [pk, ip] = max(theImpulseResponseEnvelope);

    if (options.trim)
        mask = theImpulseResponseEnvelope >= options.trim_thr * pk;
        idx  = find(mask);
        theImpulseResponseFunction = theImpulseResponseFunction(idx(1):idx(end));
        temporalSupportSeconds = temporalSupportSeconds(idx(1):idx(end));
        fprintf('Trimmed   : %d -> %d samples\n', N_up, length(h));
    end


    theImpulseResponseFunctionStruct.amplitude = theImpulseResponseFunction;
    theImpulseResponseFunctionStruct.temporalSupportSeconds = temporalSupportSeconds;

end 

function win = generateWindow(N, type, beta, alpha)
    switch lower(type)
        case 'none',     win = ones(N,1);
        case 'hann',     win = hann(N);
        case 'hamming',  win = hamming(N);
        case 'blackman', win = blackman(N);
        case 'kaiser',   win = kaiser(N, beta);
        case 'tukey',    win = tukeywin(N, alpha);
        otherwise,       error('Unknown window ''%s''.', type);
    end
end

