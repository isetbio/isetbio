function fitDelayHighPassLowPassCascadeModel()
% =========================================================
%  FIT v2: delay × HP^n_hp × LP^n_lp
%
%  Key improvement over v1:
%    - n_hp estimated from low-frequency phase offset (+90 deg per order)
%    - tau estimated from passband phase slope AFTER subtracting
%      HP and LP phase contributions
%    - tau grid search uses phase-weighted cost, not magnitude-dominated cost
%
%  Model:
%    H(jw) = K * exp(-j*w*tau) * HP(jw)^n_hp * LP(jw)^n_lp
%    HP(jw) = (j*w*tau_hp) / (1 + j*w*tau_hp)
%    LP(jw) = 1 / (1 + j*w*tau_lp)
% =========================================================

%% ---- DATA -------------------------------------------------------
frequency = (0 : 0.5 : 119.5)';
omega     = 2*pi*frequency;

transferFunction = 1e4 * [
  -0.0490+0.0000i;  0.0497-0.0058i;  0.0447-0.0020i;  0.0437+0.0060i;
   0.0452+0.0083i;  0.0481+0.0137i;  0.0519+0.0152i;  0.0551+0.0224i;
   0.0617+0.0195i;  0.0664+0.0227i;  0.0723+0.0234i;  0.0784+0.0240i;
   0.0817+0.0205i;  0.0853+0.0200i;  0.0899+0.0163i;  0.0971+0.0127i;
   0.1040+0.0085i;  0.1108+0.0036i;  0.1202-0.0102i;  0.1288-0.0170i;
   0.1396-0.0301i;  0.1458-0.0496i;  0.1527-0.0723i;  0.1514-0.1017i;
   0.1506-0.1282i;  0.1523-0.1517i;  0.2015-0.1133i;  0.1174-0.2189i;
   0.0996-0.2491i;  0.0574-0.2848i;  0.0379-0.3081i;  0.0038-0.3320i;
  -0.0599-0.3509i; -0.1287-0.3531i; -0.1999-0.3415i; -0.2523-0.3392i;
  -0.3293-0.3078i; -0.3918-0.2722i; -0.4391-0.2191i; -0.4912-0.1433i;
  -0.5429-0.0764i; -0.5872+0.0023i; -0.6004+0.0892i; -0.6121+0.1834i;
  -0.5555+0.3593i; -0.5333+0.4333i; -0.4670+0.5429i; -0.3996+0.6269i;
  -0.3472+0.7017i; -0.1725+0.8014i; -0.1167+0.8329i;  0.0628+0.8577i;
  -0.8149+0.2391i;  0.4075+0.8171i;  0.4846+0.8176i;  0.7227+0.6787i;
   0.7977+0.6025i;  0.9076+0.4610i;  1.0293+0.1963i;  1.0513+0.0859i;
  -0.4829+0.9663i;  1.0891-0.2627i;  1.0730-0.4265i;  0.9967-0.6175i;
   0.9232-0.7714i;  0.7570-0.9683i;  0.6303-1.0951i;  0.4416-1.2243i;
   0.1339-1.3305i; -0.0558-1.3698i; -0.1561-1.3893i; -0.3560-1.3780i;
  -0.6221-1.3617i; -0.9017-1.2131i; -1.0453-1.1334i; -1.2943-0.9149i;
  -1.4975-0.5809i; -1.6049-0.3164i; -1.6498-0.2644i; -1.6988+0.1468i;
  -1.7130+0.4583i; -1.6291+0.6258i; -1.4655+1.0280i; -1.4716+1.0333i;
  -1.0454+1.4987i; -0.8797+1.5951i; -0.4591+1.7301i; -0.1529+1.8034i;
   0.1019+1.8435i;  0.5501+1.7803i;  0.7346+1.7161i;  1.0333+1.5806i;
   1.3235+1.4109i;  1.4014+1.2995i;  1.6931+0.9327i;  1.7913-1.1692i;
   1.9254+0.3702i;  1.9961+0.2418i;  1.5118-1.4130i;  2.0347-0.4711i;
   1.9207-0.8154i;  1.7424-1.3113i;  1.5173-1.5132i;  1.1033-1.8273i;
   0.9202-1.9581i;  0.6384-2.0593i;  0.2742-2.1462i;  0.7002-1.9948i;
  -0.0221-2.1295i; -1.0196-1.8759i; -1.0706-1.8208i; -1.3603-1.5996i;
  -1.6288-1.2589i; -1.8944-0.9981i; -0.0045-2.0138i; -2.0948-0.1452i;
  -2.0481-0.0021i; -2.0594+0.2836i; -1.9410+0.6208i; -1.7435+0.9916i;
  -1.6130+1.2385i; -1.4946+1.3178i; -1.0532+1.6393i; -0.6370-1.8337i;
  -0.4325+1.8779i; -0.1749+1.9584i;  0.1583+1.9008i;  0.5900+1.7511i;
   0.7823+1.7186i;  0.9268+1.5845i;  1.2520+1.2887i;  1.4211+1.1097i;
   1.6179+0.7253i;  1.6796+0.4916i;  1.7961-0.1104i;  1.6933-0.2182i;
   1.6270-0.3805i;  1.5924-0.4260i;  1.2539-1.1607i;  1.1635-1.1473i;
   1.1064-1.1618i;  0.8950-1.3369i;  0.5963-1.4753i; -0.0487-1.5659i;
   0.0642-1.5112i; -0.2999-1.4947i; -0.4398-1.4277i; -0.7821-1.2092i;
  -0.8512-1.1437i; -1.1312-0.8705i; -1.1518-0.7717i; -1.1890-0.6341i;
  -1.2902-0.3294i; -1.3503-0.3013i; -1.3334+0.3117i; -1.2371+0.2831i;
  -1.1397+0.4749i; -1.0762+0.6815i; -0.9758+0.7641i; -0.7143+0.9750i;
  -0.6271+0.9837i; -0.2083+1.0904i; -0.3000+1.1108i;  0.1296+1.1232i;
   0.1170+1.0646i;  0.2518+1.0134i;  0.4253+0.9720i;  0.7931+0.6820i;
   0.8657+0.5207i;  0.7973+0.5249i;  0.8319+0.4004i;  0.9309+0.1987i;
   0.9022+0.1333i;  0.8891-0.1130i;  0.7636-0.3870i;  0.7064-0.4163i;
   0.5737-0.5970i;  0.5391-0.6371i;  0.5163-0.6019i;  0.4379-0.6252i;
   0.2975-0.6872i;  0.1741-0.7086i; -0.1467-0.7028i; -0.3025-0.6377i;
  -0.3736-0.5761i; -0.2943-0.5892i; -0.3408-0.5429i;  0.1725-0.5893i;
  -0.4633-0.3993i;  0.3491-0.4600i;  0.3621-0.4072i; -0.5764-0.0010i;
  -0.5392+0.1694i; -0.5155+0.1772i; -0.4419+0.3038i; -0.4071+0.3419i;
  -0.3614+0.3683i; -0.3725+0.3245i; -0.2944+0.3846i; -0.2375+0.4125i;
  -0.0037+0.4755i;  0.1012+0.4369i;  0.1532+0.4091i;  0.0876+0.4117i;
   0.2536+0.3307i;  0.2024+0.3376i;  0.2817+0.2682i;  0.2954+0.2320i;
   0.2961+0.2050i;  0.3365+0.1045i;  0.3246+0.0798i;  0.3213-0.0090i;
   0.2976-0.1144i;  0.2867-0.1185i;  0.2325+0.1349i;  0.2115-0.1992i;
   0.2128-0.1703i;  0.2412+0.0598i;  0.1406-0.2203i;  0.1152-0.2245i;
   0.0447+0.2483i;  0.0265-0.2352i; -0.0117-0.2273i; -0.1155-0.1827i;
  -0.1359-0.1640i; -0.1687-0.1271i; -0.1773-0.1010i; -0.1518-0.1117i;
  -0.1868-0.0364i; -0.1714-0.0547i; -0.1663-0.0426i; -0.1678+0.0027i;
  -0.1577+0.0117i; -0.1450+0.0572i; -0.1376+0.0523i; -0.1114+0.0938i;
  -0.0659+0.1207i; -0.0644+0.1145i; -0.0228+0.1280i; -0.0189+0.1214i];

N = numel(omega);

%% ---- MODEL FUNCTIONS -------------------------------------------
function Hm = model_unnorm(omega, tau, tau_hp, n_hp, tau_lp, n_lp)
    jw_hp = 1j * omega * tau_hp;
    jw_lp = 1j * omega * tau_lp;
    HP    = (jw_hp ./ (1 + jw_hp)) .^ n_hp;
    LP    = 1 ./ (1 + jw_lp) .^ n_lp;
    Hm    = exp(-1j * omega * tau) .* HP .* LP;
end

function K = opt_gain(Hm, TF, W)
    K = (conj(Hm)' * (W .* TF)) / (conj(Hm)' * (W .* Hm));
end

function r = residual(p, omega, TF, W)
    Hm  = model_unnorm(omega, p(1), p(2), p(3), p(4), p(5));
    K   = opt_gain(Hm, TF, W);
    err = TF - K * Hm;
    r   = [real(err) .* sqrt(W);  imag(err) .* sqrt(W)];
end

%% ---- STEP 1: n_hp FROM LOW-FREQUENCY PHASE OFFSET --------------
% At omega -> 0: HP contributes +n_hp*90 deg, LP and delay contribute ~0.
ph_data = unwrap(angle(transferFunction));
idx_lf  = frequency >= 0.5 & frequency <= 3;
ph_lf_deg = mean(ph_data(idx_lf)) * 180/pi;
n_hp_est  = max(1, min(6, round(ph_lf_deg / 90)));

fprintf('=== Phase-based initialization ===\n');
fprintf('  Low-freq phase ~ %.1f deg  =>  n_hp_est = %d\n', ph_lf_deg, n_hp_est);

%% ---- STEP 2: CORNER ESTIMATES FROM MAGNITUDE -------------------
mag_data   = abs(transferFunction);
idx_pass   = frequency >= 5 & frequency <= 110;
mag_peak   = max(mag_data(idx_pass));
thresh_3dB = mag_peak / sqrt(2);

idx_rise   = find(mag_data(2:end) >= thresh_3dB, 1, 'first') + 1;
idx_fall   = find(mag_data >= thresh_3dB, 1, 'last');
if isempty(idx_rise), idx_rise = 2; end
if isempty(idx_fall), idx_fall = N-1; end

f_hp_est   = max(frequency(idx_rise), 0.5);
f_lp_est   = min(frequency(idx_fall), 119.0);
tau_hp_est = 1 / (2*pi*f_hp_est);
tau_lp_est = 1 / (2*pi*f_lp_est);

fprintf('  f_hp_est = %.2f Hz,  f_lp_est = %.2f Hz\n', f_hp_est, f_lp_est);

%% ---- STEP 3: tau FROM RESIDUAL PHASE SLOPE ---------------------
% Strip HP and LP phase, regress remainder linearly vs omega (passband only).
idx_pb  = frequency >= f_hp_est & frequency <= f_lp_est;
om_pb   = omega(idx_pb);
ph_hp_k = n_hp_est * angle((1j*om_pb*tau_hp_est) ./ (1 + 1j*om_pb*tau_hp_est));
ph_lp_k = 1        * angle(1 ./ (1 + 1j*om_pb*tau_lp_est));
ph_res  = ph_data(idx_pb) - ph_hp_k - ph_lp_k;
coeffs  = [om_pb, ones(sum(idx_pb),1)] \ ph_res;
tau_slope = max(0, min(-coeffs(1), 0.08));

fprintf('  Phase-slope tau = %.4f ms\n', tau_slope*1e3);

%% ---- STEP 4: PHASE-WEIGHTED TAU GRID SEARCH -------------------
W_phase    = 1 ./ (mag_data.^2 + eps);
W_phase(1) = 0;
W_phase    = W_phase / sum(W_phase) * N;

tau_grid = linspace(0, 0.08, 8000);
cost_grid = arrayfun(@(tk) ...
    sum(abs(transferFunction - opt_gain( ...
        model_unnorm(omega,tk,tau_hp_est,n_hp_est,tau_lp_est,1), ...
        transferFunction, W_phase) .* ...
        model_unnorm(omega,tk,tau_hp_est,n_hp_est,tau_lp_est,1)).^2 .* W_phase), ...
    tau_grid);

[~, ibest] = min(cost_grid);
tau0 = tau_grid(ibest);
fprintf('  Grid-search tau = %.4f ms\n\n', tau0*1e3);

figure('Name','Tau grid','Position',[100 100 700 300]);
plot(tau_grid*1e3, cost_grid/max(cost_grid),'b-','LineWidth',1.5); hold on;
xline(tau0*1e3,'--r',sprintf('grid: %.3f ms',tau0*1e3));
xline(tau_slope*1e3,'--g',sprintf('slope: %.3f ms',tau_slope*1e3));
xlabel('\tau (ms)'); ylabel('Norm. cost'); title('Phase-weighted delay grid'); grid on;
legend('Cost','Grid min','Phase slope','Location','best');
saveas(gcf,'tau_grid_v2.png');

%% ---- STEP 5: FULL IRLS FIT -------------------------------------
p0 = [0.02,        tau_hp_est,   n_hp_est,  tau_lp_est,  1.0];
lb = [0.0,        1/(2*pi*119), 0.5,        1/(2*pi*119), 0.5];
ub = [0.06,        1/(2*pi*0.5), 2,        1/(2*pi*0.5), 2];

opts = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...
    'MaxFunEvals',20000,'MaxIterations',3000, ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
    'Display','iter');

W      = ones(N,1);
W(1)   = 0;

for iter = 1:8
    fprintf('=== IRLS %d/8 ===\n', iter);
    p_fit = lsqnonlin(@(p) residual(p, omega, transferFunction, W), p0, lb, ub, opts);
    Hm    = model_unnorm(omega, p_fit(1),p_fit(2),p_fit(3),p_fit(4),p_fit(5));
    K     = opt_gain(Hm, transferFunction, W);
    err   = transferFunction - K*Hm;
    u     = abs(err) / (4.685 * median(abs(err(W>0))) + eps);
    W     = (u < 1) .* (1 - u.^2).^2;
    W(1)  = 0;
    W     = W / (sum(W)+eps) * N;
    p0    = p_fit;
end

%% ---- REPORT ----------------------------------------------------
tau=p_fit(1); tau_hp=p_fit(2); n_hp=p_fit(3); tau_lp=p_fit(4); n_lp=p_fit(5);
f_hp = 1/(2*pi*tau_hp);  f_lp = 1/(2*pi*tau_lp);
Hm    = model_unnorm(omega, tau, tau_hp, n_hp, tau_lp, n_lp);
K     = opt_gain(Hm, transferFunction, W);
H_fit = K * Hm;
resid = transferFunction - H_fit;
R2    = 1 - sum(abs(resid).^2) / sum(abs(transferFunction-mean(transferFunction)).^2);
wRMSE = sqrt(sum(W.*abs(resid).^2)/sum(W));
outlier_mask = W < 0.3;

fprintf('\n=== FITTED PARAMETERS ===\n');
fprintf('  tau       = %.4f ms\n',  tau*1e3);
fprintf('  f_hp      = %.3f Hz  (n_hp = %.3f)\n', f_hp, n_hp);
fprintf('  f_lp      = %.3f Hz  (n_lp = %.3f)\n', f_lp, n_lp);
fprintf('  |K|       = %.4e,  angle(K) = %.2f deg\n', abs(K), angle(K)*180/pi);
fprintf('  R2        = %.6f\n', R2);
fprintf('  wRMSE     = %.4e\n', wRMSE);
fprintf('  Outliers  = %d / %d\n', sum(outlier_mask), N);

%% ---- PLOTS -----------------------------------------------------
f_d  = linspace(0.1,120,3000)';
w_d  = 2*pi*f_d;
Hm_d = K * model_unnorm(w_d, tau, tau_hp, n_hp, tau_lp, n_lp);

fig = figure('Name','Fit v2','Position',[80 80 1200 850]);

subplot(2,3,1);
semilogy(frequency,abs(transferFunction),'o','MarkerSize',3,'Color',[.5 .5 .5],'DisplayName','Data'); hold on;
semilogy(f_d,abs(Hm_d),'r-','LineWidth',2,'DisplayName','Fit');
semilogy(frequency(outlier_mask),abs(transferFunction(outlier_mask)),'x','MarkerSize',8,'Color',[.85 .5 0],'LineWidth',1.5,'DisplayName','Outlier');
xline(f_hp,'--b',sprintf('f_{hp}=%.1f Hz',f_hp),'LabelVerticalAlignment','bottom');
xline(f_lp,'--m',sprintf('f_{lp}=%.1f Hz',f_lp),'LabelVerticalAlignment','bottom');
xlabel('Frequency (Hz)'); ylabel('|H|'); title('Magnitude'); legend('Location','best'); grid on;

subplot(2,3,2);
ph_uw = unwrap(angle(transferFunction));
plot(frequency,ph_uw*180/pi,'o','MarkerSize',3,'Color',[.5 .5 .5],'DisplayName','Data'); hold on;
plot(f_d,unwrap(angle(Hm_d))*180/pi,'r-','LineWidth',2,'DisplayName','Fit');
plot(frequency(outlier_mask),ph_uw(outlier_mask)*180/pi,'x','MarkerSize',8,'Color',[.85 .5 0],'LineWidth',1.5,'DisplayName','Outlier');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); title('Phase (unwrapped)'); legend('Location','best'); grid on;

subplot(2,3,3);
plot(real(transferFunction),imag(transferFunction),'o','MarkerSize',3,'Color',[.5 .5 .5],'DisplayName','Data'); hold on;
plot(real(Hm_d),imag(Hm_d),'r-','LineWidth',2,'DisplayName','Fit');
plot(real(transferFunction(outlier_mask)),imag(transferFunction(outlier_mask)),'x','MarkerSize',8,'Color',[.85 .5 0],'LineWidth',1.5,'DisplayName','Outlier');
xlabel('Real'); ylabel('Imag'); title('Nyquist'); legend('Location','best'); grid on; axis equal;

subplot(2,3,4);
stem(frequency,abs(resid),'filled','MarkerSize',3,'Color',[.3 .3 .7]); hold on;
stem(frequency(outlier_mask),abs(resid(outlier_mask)),'filled','MarkerSize',5,'Color',[.85 .5 0]);
xlabel('Frequency (Hz)'); ylabel('|residual|'); title('Residual magnitude'); grid on;

subplot(2,3,5);
stem(frequency,W,'filled','MarkerSize',3,'Color',[.2 .6 .4]);
yline(0.3,'--r','Outlier threshold');
xlabel('Frequency (Hz)'); ylabel('Weight'); title('IRLS weights'); grid on; ylim([0,max(W)*1.1]);

subplot(2,3,6);
plot(frequency,real(transferFunction),'o','MarkerSize',2,'Color',[.4 .6 .9],'DisplayName','Re data'); hold on;
plot(f_d,real(Hm_d),'-','Color',[.1 .3 .8],'LineWidth',1.5,'DisplayName','Re fit');
plot(frequency,imag(transferFunction),'o','MarkerSize',2,'Color',[.9 .6 .4],'DisplayName','Im data');
plot(f_d,imag(Hm_d),'-','Color',[.8 .3 .1],'LineWidth',1.5,'DisplayName','Im fit');
xlabel('Frequency (Hz)'); ylabel('H'); title('Real & Imaginary'); legend('Location','best'); grid on;

sgtitle(sprintf('HP×LP×Delay v2  |  \\tau=%.3f ms  f_{hp}=%.1f Hz (n_{hp}=%.2f)  f_{lp}=%.1f Hz (n_{lp}=%.2f)  R^2=%.5f', ...
    tau*1e3, f_hp, n_hp, f_lp, n_lp, R2),'FontSize',10);

saveas(fig,'fit_hp_lp_delay_v2.png');
fprintf('\nSaved: tau_grid_v2.png,  fit_hp_lp_delay_v2.png\n');
end