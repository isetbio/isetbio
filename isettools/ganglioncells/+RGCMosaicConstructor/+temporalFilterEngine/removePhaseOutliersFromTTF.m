%
% RGCMosaicConstructor.temporalFilterEngine.removePhaseOutliersFromTTF(temporalFrequecySupport, theTTF)
%

function theCleanedTTF = removePhaseOutliersFromTTF(temporalFrequecySupport, theTTF)
% =========================================================
%  PHASE SPECTRUM OUTLIER DETECTION
%  No toolbox dependencies (base MATLAB only)
%
%  Expects in workspace:
%    temporalFrequecySupport        - vector of frequencies (Hz), any spacing
%    transferFunction - complex vector, same length as temporalFrequecySupport
%
%  Methods:
%    1. Velocity  (1st finite difference) - abrupt jumps
%    2. Acceleration (2nd finite difference) - curvature spikes
%    3. Robust Whittaker smoother residual 
%    4. MAD on 2nd differences (independent check on method 2)
%
%  The Whittaker smoother solves:
%    min  sum(w*(y - z)^2) + lambda * sum((D2*z)^2)
%  where D2 is the 2nd-order difference matrix. This is a
%  pure linear algebra problem (sparse backslash), requiring
%  no toolboxes.
%
%  A point is flagged as outlier if >= min_votes methods agree.
%
%  Outputs:
%    outlier_mask           - logical N x 1, true = outlier
%    theCleanedTTF - TF with outlier phases interpolated
% =========================================================

beVerbose = false;

% ---- TUNING PARAMETERS ----------------------------------------
k_vel     = 5;       % Method 1: MAD multiplier for 1st-diff threshold
k_acc     = 5;       % Method 2: MAD multiplier for 2nd-diff threshold
k_spl     = 4;       % Method 3: MAD multiplier for smoother residual
k_mad2    = 5;       % Method 4: MAD multiplier for 2nd-diff (independent)
lambda    = 1e4;     % Whittaker smoothing strength (larger = smoother)
min_votes = 2;       % Min methods that must agree to call outlier

% UNWRAP PHASE
ph = unwrap(angle(theTTF(:)));
f  = temporalFrequecySupport(:);
N  = numel(ph);

% ---- METHOD 1: VELOCITY (1st finite difference) ----------------
d1        = diff(ph);                   % (N-1) x 1
med_d1    = median(d1);
mad_d1    = 1.4826 * median(abs(d1 - med_d1));
thr1      = k_vel * mad_d1;
bad_steps = abs(d1 - med_d1) > thr1;

flag1 = false(N, 1);
flag1(1:N-1) = flag1(1:N-1) | bad_steps;
flag1(2:N)   = flag1(2:N)   | bad_steps;

% ---- METHOD 2: ACCELERATION (2nd finite difference) -----------
d2      = diff(ph, 2);                  % (N-2) x 1
med_d2  = median(d2);
mad_d2  = 1.4826 * median(abs(d2 - med_d2));
thr2    = k_acc * mad_d2;
bad_acc = abs(d2 - med_d2) > thr2;

flag2 = false(N, 1);
flag2(2:N-1) = flag2(2:N-1) | bad_acc;

% Method 3 Iteratively re-weighted Whittaker (bisquare weights)
W_spl = ones(N, 1);

for k = 1:6
    ph_smooth = whittaker(ph, W_spl, lambda, N);
    res_spl   = ph - ph_smooth;
    mad_spl   = 1.4826 * median(abs(res_spl - median(res_spl)));
    u         = res_spl / (4.685 * mad_spl + eps);
    W_spl     = (abs(u) < 1) .* (1 - u.^2).^2;
    W_spl     = max(W_spl, 1e-6);
end
thr3  = k_spl * mad_spl;
flag3 = abs(res_spl - median(res_spl)) > thr3;

% METHOD 4: MAD ON 2ND DIFFERENCES (independent) ----------
d2b     = diff(ph, 2);
mad_d2b = 1.4826 * median(abs(d2b));
thr4    = k_mad2 * mad_d2b;
flag4   = false(N, 1);
flag4(2:N-1) = abs(d2b) > thr4;

% Vote on the 4 methods to determine the outliers
votes        = int32(flag1) + int32(flag2) + int32(flag3) + int32(flag4);
outlier_mask = votes >= min_votes;

% Interpolate phase at outlier points
ph_clean  = ph;
clean_mask = ~outlier_mask;
if any(outlier_mask)
    ph_clean(outlier_mask) = interp1( ...
        f(clean_mask), ph(clean_mask), ...
        f(outlier_mask), 'pchip');
end

% The cleaned TTF
theCleanedTTF = reshape(abs(theTTF(:)) .* exp(1j * ph_clean), size(theTTF));


if (beVerbose)
    outlier_idx  = find(outlier_mask);
    fprintf('=== PHASE OUTLIER DETECTION ===\n');
    fprintf('  Thresholds: vel=%.3f rad, acc=%.3f rad, smooth=%.3f rad, mad2=%.3f rad\n', ...
        thr1, thr2, thr3, thr4);
    fprintf('  Total outliers (%d+ / 4 votes): %d / %d\n\n', min_votes, numel(outlier_idx), N);
    fprintf('  Index   Freq(Hz)   Phase(deg)   Votes\n');
    fprintf('  ------  --------   ----------   -----\n');
    for i = 1:numel(outlier_idx)
        ii = outlier_idx(i);
        fprintf('  %4d    %6.1f     %9.2f    %d\n', ii, f(ii), ph(ii)*180/pi, votes(ii));
    end
    fprintf('\n  Clean points: %d / %d\n', N - numel(outlier_idx), N);
end




end


function z = whittaker(ph, W_vec, lambda, N)

    % WHITTAKER SMOOTHER 
    % Builds a sparse 2nd-order difference matrix D2 and solves
    % the penalized least squares problem iteratively for robustness.
    %
    % D2 is (N-2) x N such that D2*z = diff(z, 2).
    e  = ones(N, 1);
    D1 = spdiags([-e, e],  0:1, N-1, N);
    D2 = spdiags([ e, -2*e, e], 0:2, N-2, N);   % 2nd-order difference matrix


    % Solve: (diag(W) + lambda * D2' * D2) * z = W .* ph
    W_diag = spdiags(W_vec, 0, numel(ph), numel(ph));
    A      = W_diag + lambda * (D2' * D2);
    z      = A \ (W_vec .* ph);
end
