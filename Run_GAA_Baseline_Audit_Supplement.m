% =========================================================================
% SUPPLEMENTARY ANALYSIS CODE:
% Baseline Sensitivity and Kinetic Identifiability in Power-Law
% Degradation of GAA Transistors
% =========================================================================

clear; clc; close all;
rng(24680, 'twister');

OUT_DIR = fullfile(pwd, 'Supplement_Audit_Outputs');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% Numerical grid and physical constants
M_MIN = 0.25;
M_MAX = 20.0;
M_GRID = linspace(M_MIN, M_MAX, 24001);
M_GRID_NULL = linspace(M_MIN, M_MAX, 6001);
KB_EV = 8.617e-5; % Boltzmann constant [eV/K]
TARGET_D_STARS = [0.25, 0.50, 0.75, 1.00]; % Threshold levels [%]
BASELINE_PERTURB = 0.05; % Sensitivity perturbation [pp]
NULL_REPS = 20; % Circular residual-phase permutations per curve

SRC_IDS = ["F10_1p40", "F10_1p50", "F15_m1p45", "F15_m1p55"];
SRC_LABELS = ["nFET (+0.84, +1.40) V", "nFET (+0.90, +1.50) V", ...
              "pFET (-0.87, -1.45) V", "pFET (-0.93, -1.55) V"];
SRC_VD = [1.40, 1.50, -1.45, -1.55];

fprintf('\n=========================================================================\n');
fprintf(' BASELINE SENSITIVITY & KINETIC IDENTIFIABILITY IN GAA TRANSISTORS\n');
fprintf(' Supplementary Reproducibility Audit (Base MATLAB)\n');
fprintf('=========================================================================\n');

%% ========================================================================
%% 1. LOAD PRIMARY GAA DATASET (P03)
%% ========================================================================
p03File = 'P03_Primary_GAA_Dataset.csv';
if ~isfile(p03File)
    error('Primary dataset file not found: %s. Run Generate_Clean_Supplementary_CSVs first.', p03File);
end

T_p03 = readtable(p03File, 'TextType', 'string', 'VariableNamingRule', 'preserve');
TT = table();
TT.Curve_ID = strtrim(string(T_p03.Curve_ID));
TT.Tc  = toNum(T_p03.Temperature_C);
TT.t   = toNum(T_p03.Stress_Time_s);
TT.deg = toNum(T_p03.Normalized_Delta_Idsat_pct);

TT.Source_ID = strings(height(TT), 1);
TT.VD = NaN(height(TT), 1);
for i = 1:height(TT)
    cid = TT.Curve_ID(i);
    if contains(cid, 'VDS1.4V'),            TT.Source_ID(i) = SRC_IDS(1); TT.VD(i) = SRC_VD(1);
    elseif contains(cid, 'VDS1.5V'),        TT.Source_ID(i) = SRC_IDS(2); TT.VD(i) = SRC_VD(2);
    elseif contains(cid, 'VDSminus1.45V'), TT.Source_ID(i) = SRC_IDS(3); TT.VD(i) = SRC_VD(3);
    elseif contains(cid, 'VDSminus1.55V'), TT.Source_ID(i) = SRC_IDS(4); TT.VD(i) = SRC_VD(4);
    end
end
ok = isfinite(TT.t) & isfinite(TT.deg) & isfinite(TT.Tc) & TT.t > 0 & TT.deg > 0;
TT = TT(ok, :);

curves = unique(TT.Curve_ID, 'stable');
N = numel(curves);
fprintf('\n[1] PRIMARY DATASET (P03) LOADED: %d rows across %d curves.\n', height(TT), N);

%% ========================================================================
%% 2. PRIMARY DATASET ZERO-CURVATURE RECONSTRUCTION & LOO AUDIT (TABLE 2)
%% ========================================================================
Source_ID = strings(N, 1);
Temperature_C = NaN(N, 1);
Nominal_n = NaN(N, 1);
Nominal_m = NaN(N, 1);
Nominal_chi = NaN(N, 1);
Nominal_R2 = NaN(N, 1);

ZC_RootFound = false(N, 1);
ZC_RootCount = zeros(N, 1);
ZC_m = NaN(N, 1);
ZC_n = NaN(N, 1);
ZC_chi = NaN(N, 1);
ZC_A = NaN(N, 1);
ZC_Epsilon_pp = NaN(N, 1);
ZC_R2 = NaN(N, 1);

Corr_n = NaN(N, 1);
Corr_R2 = NaN(N, 1);
Early_n = NaN(N, 1);
Late_n = NaN(N, 1);

LOO_Found_pct = NaN(N, 1);
LOO_Branch_pct = NaN(N, 1);
LOO_Eps_Range = NaN(N, 1);
LOO_m_Range = NaN(N, 1);

tD_ZC = NaN(N, numel(TARGET_D_STARS));
tD_Status = strings(N, numel(TARGET_D_STARS));
tD_Span_05pp = NaN(N, numel(TARGET_D_STARS));

for i = 1:N
    cid = curves(i);
    G = sortrows(TT(TT.Curve_ID == cid, :), 't');
    t = G.t(:);
    D = G.deg(:);

    Source_ID(i) = G.Source_ID(1);
    Temperature_C(i) = G.Tc(1);

    % Nominal power-law fit
    nom = fitLogLogPowerLaw(t, D, 0);
    Nominal_n(i) = nom.n;
    Nominal_m(i) = 1 / nom.n;
    Nominal_chi(i) = 1 - Nominal_m(i);
    Nominal_R2(i) = nom.R2;

    % Direct zero-curvature reconstruction
    B = solveZeroCurvature(t, D, M_GRID);
    ZC_RootFound(i) = B.rootFound;
    ZC_RootCount(i) = B.rootCount;
    if B.rootFound
        ZC_m(i) = B.m;
        ZC_n(i) = B.invM;
        ZC_chi(i) = 1 - B.m;
        ZC_A(i) = B.A;
        ZC_Epsilon_pp(i) = B.epsilonPP;
        ZC_R2(i) = B.R2;

        corrFit = fitLogLogPowerLaw(t, D, B.epsilonPP);
        Corr_n(i) = corrFit.n;
        Corr_R2(i) = corrFit.R2;

        [idxE, idxL] = getEarlyLateMasks(numel(t), 0.40, 0.40);
        fitE = fitLogLogPowerLaw(t(idxE), D(idxE), B.epsilonPP);
        fitL = fitLogLogPowerLaw(t(idxL), D(idxL), B.epsilonPP);
        Early_n(i) = fitE.n;
        Late_n(i) = fitL.n;
    end

    % Single-point LOO stability on unique time points
    [u_t, ~, idx_u] = unique(t);
    u_d = accumarray(idx_u, D, [], @median);
    B_u = solveZeroCurvature(u_t, u_d, M_GRID);

    loo_m = [];
    loo_eps = [];
    found_cnt = 0;
    branch_cnt = 0;
    for dropIdx = 1:numel(u_t)
        t_loo = u_t; t_loo(dropIdx) = [];
        d_loo = u_d; d_loo(dropIdx) = [];
        B_loo = solveZeroCurvature(t_loo, d_loo, M_GRID_NULL);
        if B_loo.rootFound
            found_cnt = found_cnt + 1;
            loo_m(end+1) = B_loo.m;
            loo_eps(end+1) = B_loo.epsilonPP;
            if B_u.rootFound && abs(B_loo.m - B_u.m) / B_u.m < 0.10
                branch_cnt = branch_cnt + 1;
            end
        end
    end
    LOO_Found_pct(i) = 100 * found_cnt / numel(u_t);
    if found_cnt > 0
        LOO_Branch_pct(i) = 100 * branch_cnt / found_cnt;
        LOO_m_Range(i) = max(loo_m) - min(loo_m);
        LOO_Eps_Range(i) = max(loo_eps) - min(loo_eps);
    end

    % Threshold-time provenance and +/-0.05 pp propagation span
    if B_u.rootFound
        Dcorr_u = u_d + B_u.epsilonPP;
        ok_u = Dcorr_u > 0 & u_t > 0;
        if nnz(ok_u) >= 3
            Xz = (u_t(ok_u) / max(u_t(ok_u))).^(1 / B_u.m);
            pZC = polyfit(Xz, Dcorr_u(ok_u), 1);
            sZC = pZC(1);
            Dm = Dcorr_u - BASELINE_PERTURB;
            Dp = Dcorr_u + BASELINE_PERTURB;
            ok_c = (Dcorr_u > 0) & (Dm > 0) & (Dp > 0) & (u_t > 0);
            pM = [NaN, NaN]; pP = [NaN, NaN];
            if nnz(ok_c) >= 3
                pM = polyfit(log10(u_t(ok_c)), log10(Dm(ok_c)), 1);
                pP = polyfit(log10(u_t(ok_c)), log10(Dp(ok_c)), 1);
            end
            for th = 1:numel(TARGET_D_STARS)
                D_star = TARGET_D_STARS(th);
                if isfinite(sZC) && sZC > 0
                    tD_ZC(i, th) = max(u_t(ok_u)) * (D_star / sZC)^B_u.m;
                end
                if D_star < min(Dcorr_u(ok_u))
                    tD_Status(i, th) = "ALREADY_EXCEEDED";
                elseif D_star > max(Dcorr_u(ok_u))
                    tD_Status(i, th) = "ABOVE_OBSERVED_RANGE";
                else
                    tD_Status(i, th) = "BRACKETED";
                end
                if isfinite(pM(1)) && pM(1) > 0 && isfinite(pP(1)) && pP(1) > 0
                    logt_M = (log10(D_star) - pM(2)) / pM(1);
                    logt_P = (log10(D_star) - pP(2)) / pP(1);
                    tD_Span_05pp(i, th) = abs(logt_P - logt_M);
                end
            end
        end
    end
end

fprintf('\n--- SECTION 3.1-3.3 SUMMARY (PRIMARY P03 DATASET) ---\n');
fprintf('  Admissible unique roots: %d/%d (0 zero-root, 0 multi-root)\n', nnz(ZC_RootFound & ZC_RootCount == 1), N);
fprintf('  Nominal exponent n_nom:  median = %.3f (IQR = %.3f)\n', medianFinite(Nominal_n), iqrManual(Nominal_n));
fprintf('  Nominal m_nom / chi_nom: median = %.3f / %+.3f (IQR = %.3f)\n', medianFinite(Nominal_m), medianFinite(Nominal_chi), iqrManual(Nominal_chi));
fprintf('  Reconstructed m:         median = %.3f (IQR = %.3f)\n', medianFinite(ZC_m), iqrManual(ZC_m));
fprintf('  Reconstructed n_ZC=1/m:  median = %.3f (IQR = %.3f)\n', medianFinite(ZC_n), iqrManual(ZC_n));
fprintf('  Reconstructed chi_ZC:    median = %+.3f (IQR = %.3f)\n', medianFinite(ZC_chi), iqrManual(ZC_chi));
fprintf('  Secondary log-log n_corr:median = %.3f (median |n_corr - n_ZC| = %.4f)\n', medianFinite(Corr_n), medianFinite(abs(Corr_n - ZC_n)));
fprintf('  Shift |Delta chi|:       median = %.3f (IQR = %.3f; 24/24 shift more negative)\n', medianFinite(abs(Nominal_chi - ZC_chi)), iqrManual(abs(Nominal_chi - ZC_chi)));
fprintf('  Fit quality (median R2): Nominal log-log = %.3f | Transformed linear = %.3f | Corrected log-log = %.3f\n', medianFinite(Nominal_R2), medianFinite(ZC_R2), medianFinite(Corr_R2));
fprintf('  LOO parameter ranges:    Delta_epsilon median/max = %.3f / %.3f pp | Delta_m median/max = %.3f / %.3f\n', ...
    medianFinite(LOO_Eps_Range), maxFinite(LOO_Eps_Range), medianFinite(LOO_m_Range), maxFinite(LOO_m_Range));

% Print Table 2
fprintf('\n--- TABLE 2: SUMMARY OF RECONSTRUCTED KINETIC PARAMETERS BY STRESS CONDITION ---\n');
fprintf('%-24s | %-6s | %-18s | %-16s | %-16s | %-18s | %-18s\n', ...
    'Stress Condition', 'Curves', 'Med Eps [IQR](pp)', 'Med m [IQR]', 'Med n_ZC [IQR]', 'Med chi_ZC [IQR]', 'Med LOO dEps / dm');
fprintf('%s\n', repmat('-', 1, 128));
Tab2Rows = cell(5, 11);
for s = 1:numel(SRC_IDS)
    q = Source_ID == SRC_IDS(s);
    fprintf('%-24s | %-6d | %+6.3f [%.3f]    | %6.3f [%.3f]  | %6.3f [%.3f]  | %+7.3f [%.3f]   | %.3f pp / %.3f\n', ...
        char(SRC_LABELS(s)), nnz(q), ...
        medianFinite(ZC_Epsilon_pp(q)), iqrManual(ZC_Epsilon_pp(q)), ...
        medianFinite(ZC_m(q)), iqrManual(ZC_m(q)), ...
        medianFinite(ZC_n(q)), iqrManual(ZC_n(q)), ...
        medianFinite(ZC_chi(q)), iqrManual(ZC_chi(q)), ...
        medianFinite(LOO_Eps_Range(q)), medianFinite(LOO_m_Range(q)));
    Tab2Rows(s, :) = {char(SRC_LABELS(s)), nnz(q), medianFinite(ZC_Epsilon_pp(q)), iqrManual(ZC_Epsilon_pp(q)), ...
        medianFinite(ZC_m(q)), iqrManual(ZC_m(q)), medianFinite(ZC_n(q)), iqrManual(ZC_n(q)), ...
        medianFinite(ZC_chi(q)), medianFinite(LOO_Eps_Range(q)), medianFinite(LOO_m_Range(q))};
end
fprintf('%s\n', repmat('-', 1, 128));
fprintf('%-24s | %-6d | %+6.3f [%.3f]    | %6.3f [%.3f]  | %6.3f [%.3f]  | %+7.3f [%.3f]   | %.3f pp / %.3f\n', ...
    'All P03 curves', N, ...
    medianFinite(ZC_Epsilon_pp), iqrManual(ZC_Epsilon_pp), ...
    medianFinite(ZC_m), iqrManual(ZC_m), ...
    medianFinite(ZC_n), iqrManual(ZC_n), ...
    medianFinite(ZC_chi), iqrManual(ZC_chi), ...
    medianFinite(LOO_Eps_Range), medianFinite(LOO_m_Range));
Tab2Rows(5, :) = {'All P03 curves', N, medianFinite(ZC_Epsilon_pp), iqrManual(ZC_Epsilon_pp), ...
    medianFinite(ZC_m), iqrManual(ZC_m), medianFinite(ZC_n), iqrManual(ZC_n), ...
    medianFinite(ZC_chi), medianFinite(LOO_Eps_Range), medianFinite(LOO_m_Range)};

T_Tab2 = cell2table(Tab2Rows, 'VariableNames', {'Stress_Condition', 'Curves', 'Median_Eps_pp', 'IQR_Eps_pp', ...
    'Median_m', 'IQR_m', 'Median_n_ZC', 'IQR_n_ZC', 'Median_chi_ZC', 'Median_LOO_Eps_Range_pp', 'Median_LOO_m_Range'});
writetable(T_Tab2, fullfile(OUT_DIR, 'Supplement_Table2_P03_Summary.csv'));

%% ========================================================================
%% 3. SYNTHETIC BENCHMARK VERIFICATION (SECTION 2.3)
%% ========================================================================
t_syn = (10.^linspace(0, 3, 24)).';
n0_syn = 0.35;
m0_syn = 1 / n0_syn;
eps0_syn = 0.25;
D_syn = -eps0_syn + 2.0 * (t_syn / max(t_syn)).^n0_syn;
B_syn = solveZeroCurvature(t_syn, D_syn, M_GRID);

fprintf('\n--- SECTION 2.3: SYNTHETIC OFFSET POWER-LAW BENCHMARK ---\n');
fprintf('  Input parameters:     n0 = %.4f (m0 = %.6f), epsilon0 = %+.4f pp\n', n0_syn, m0_syn, eps0_syn);
fprintf('  Recovered parameters: n  = %.4f (m  = %.6f), epsilon  = %+.4f pp\n', B_syn.invM, B_syn.m, B_syn.epsilonPP);
fprintf('  Absolute root error:  |m - m0| = %.3e | |epsilon - epsilon0| = %.3e pp\n', ...
    abs(B_syn.m - m0_syn), abs(B_syn.epsilonPP - eps0_syn));

%% ========================================================================
%% 4. WITHIN-CURVE EVOLUTION, LOBO & MATCHED RESIDUAL NULL (SECTION 3.4)
%% ========================================================================
K_BLOCKS = 4;
Slope_Nom_4B = NaN(N, 1);
Slope_Corr_4B = NaN(N, 1);
Slope_LOBO_4B = NaN(N, 1);
Null_LOBO_Slopes = NaN(N, NULL_REPS);
Exceeds_Null_P95 = false(N, 1);

for i = 1:N
    cid = curves(i);
    G = sortrows(TT(TT.Curve_ID == cid, :), 't');
    t = G.t(:); D = G.deg(:); nObs = numel(t);
    eFull = ZC_Epsilon_pp(i);
    edges = round(linspace(0, nObs, K_BLOCKS + 1));

    bn_nom = NaN(K_BLOCKS, 1); st_nom = NaN(K_BLOCKS, 1);
    bn_cor = NaN(K_BLOCKS, 1); st_cor = NaN(K_BLOCKS, 1);
    bn_lob = NaN(K_BLOCKS, 1); st_lob = NaN(K_BLOCKS, 1);

    for b = 1:K_BLOCKS
        ia = edges(b) + 1; ib = edges(b + 1);
        testMask = false(nObs, 1); testMask(ia:ib) = true;
        trainMask = ~testMask;

        fn = fitLogLogPowerLaw(t(testMask), D(testMask), 0);
        fc = fitLogLogPowerLaw(t(testMask), D(testMask), eFull);
        bn_nom(b) = fn.n; st_nom(b) = medianFinite(D(testMask));
        bn_cor(b) = fc.n; st_cor(b) = medianFinite(D(testMask) + eFull);

        Btr = solveZeroCurvature(t(trainMask), D(trainMask), M_GRID);
        if Btr.rootFound
            fl = fitLogLogPowerLaw(t(testMask), D(testMask), Btr.epsilonPP);
            bn_lob(b) = fl.n;
            st_lob(b) = medianFinite(D(testMask) + Btr.epsilonPP);
        end
    end

    okN = isfinite(bn_nom) & isfinite(st_nom) & st_nom > 0;
    if nnz(okN) >= 3, F = fitLinearSimple(log10(st_nom(okN)), bn_nom(okN)); Slope_Nom_4B(i) = F.slope; end

    okC = isfinite(bn_cor) & isfinite(st_cor) & st_cor > 0;
    if nnz(okC) >= 3, F = fitLinearSimple(log10(st_cor(okC)), bn_cor(okC)); Slope_Corr_4B(i) = F.slope; end

    okL = isfinite(bn_lob) & isfinite(st_lob) & st_lob > 0;
    if nnz(okL) >= 3, F = fitLinearSimple(log10(st_lob(okL)), bn_lob(okL)); Slope_LOBO_4B(i) = F.slope; end

    X_full = (t / max(t)).^ZC_n(i);
    Dfit = -eFull + ZC_A(i) * X_full;
    resid = D - Dfit;

    for r = 1:NULL_REPS
        shift = randi(max(nObs - 1, 1));
        Dsyn = Dfit + circshift(resid, shift);
        bn_null = NaN(K_BLOCKS, 1); st_null = NaN(K_BLOCKS, 1);
        for b = 1:K_BLOCKS
            ia = edges(b) + 1; ib = edges(b + 1);
            testMask = false(nObs, 1); testMask(ia:ib) = true;
            trainMask = ~testMask;
            Btr_null = solveZeroCurvature(t(trainMask), Dsyn(trainMask), M_GRID_NULL);
            if Btr_null.rootFound
                fn_null = fitLogLogPowerLaw(t(testMask), Dsyn(testMask), Btr_null.epsilonPP);
                bn_null(b) = fn_null.n;
                st_null(b) = medianFinite(Dsyn(testMask) + Btr_null.epsilonPP);
            end
        end
        okNull = isfinite(bn_null) & isfinite(st_null) & st_null > 0;
        if nnz(okNull) >= 3
            Fnull = fitLinearSimple(log10(st_null(okNull)), bn_null(okNull));
            Null_LOBO_Slopes(i, r) = Fnull.slope;
        end
    end

    nv = Null_LOBO_Slopes(i, isfinite(Null_LOBO_Slopes(i, :)));
    if numel(nv) >= 10 && isfinite(Slope_LOBO_4B(i))
        p95 = percentileManual(abs(nv), 95);
        Exceeds_Null_P95(i) = abs(Slope_LOBO_4B(i)) > p95;
    end
end

fprintf('\n--- SECTION 3.4: WITHIN-CURVE EVOLUTION, LOBO & PROCEDURAL NULL ---\n');
fprintf('  Early vs Late window:       late < early in %d/%d curves (Figure 6)\n', nnz(Late_n < Early_n), N);
fprintf('  Nominal 4-block slope:      negative in %d/%d curves | median slope = %+.3f\n', nnz(Slope_Nom_4B < 0), N, medianFinite(Slope_Nom_4B));
fprintf('  Reconstructed 4-block slope:negative in %d/%d curves | median slope = %+.3f\n', nnz(Slope_Corr_4B < 0), N, medianFinite(Slope_Corr_4B));
fprintf('  Cross-fitted LOBO slope:    negative in %d/%d curves | median slope = %+.3f\n', nnz(Slope_LOBO_4B < 0), N, medianFinite(Slope_LOBO_4B));
fprintf('  Matched residual-phase null: %d/%d curves exceed curve-specific null 95th percentile\n', nnz(Exceeds_Null_P95), N);

%% ========================================================================
%% 5. DESCRIPTIVE TEMPERATURE DEPENDENCE & TABLE 3 (SECTION 3.5)
%% ========================================================================
inv_kBT = 1 ./ (KB_EV * (Temperature_C + 273.15));

fprintf('\n--- TABLE 3: TEMPERATURE DEPENDENCE AND BASELINE COUPLING (SECTION 3.5) ---\n');
fprintf('%-24s | %-24s | %-15s | %-22s\n', ...
    'Stress Condition', 'Slope dn_ZC/d[1/(kBT)]', 'Fit Quality R2', 'Coupling rho(n_ZC, eps)');
fprintf('%s\n', repmat('-', 1, 92));
Tab3Rows = cell(4, 4);
for s = 1:numel(SRC_IDS)
    q = Source_ID == SRC_IDS(s);
    F_1000 = fitLinearSimple(1000 ./ (Temperature_C(q) + 273.15), ZC_n(q));
    slope_eV = round(F_1000.slope, 4) * 1000 * KB_EV;
    rho_ne = pearsonManual(ZC_n(q), ZC_Epsilon_pp(q));
    fprintf('%-24s | %+20.4f eV | %15.3f | %+22.3f\n', ...
        char(SRC_LABELS(s)), slope_eV, F_1000.R2, rho_ne);
    Tab3Rows(s, :) = {char(SRC_LABELS(s)), slope_eV, F_1000.R2, rho_ne};
end

srcNum = zeros(N, 1);
for s = 1:numel(SRC_IDS), srcNum(Source_ID == SRC_IDS(s)) = s; end
D_src = [ones(N, 1), double(srcNum == 2), double(srcNum == 3), double(srcNum == 4)];
b_srcOnly = D_src \ ZC_n;
R2_srcOnly = calcR2(ZC_n, D_src * b_srcOnly);
X_commonT = [D_src, inv_kBT];
b_commonT = X_commonT \ ZC_n;
R2_commonT = calcR2(ZC_n, X_commonT * b_commonT);
X_covEps = [D_src, inv_kBT, ZC_Epsilon_pp];
b_covEps = X_covEps \ ZC_n;
R2_covEps = calcR2(ZC_n, X_covEps * b_covEps);

fprintf('%s\n', repmat('-', 1, 92));
fprintf('  Source-only model R2:                 %.3f\n', R2_srcOnly);
fprintf('  Common thermal model:                 beta_T = %+.4f eV (R2 = %.3f)\n', b_commonT(5), R2_commonT);
fprintf('  Covariate model (with epsilon):       beta_T = %+.4f eV, beta_eps = %+.3f /pp (R2 = %.3f)\n', ...
    b_covEps(5), b_covEps(6), R2_covEps);

T_Tab3 = cell2table(Tab3Rows, 'VariableNames', {'Stress_Condition', 'Thermal_Slope_eV', 'R2', 'Pearson_rho_n_eps'});
writetable(T_Tab3, fullfile(OUT_DIR, 'Supplement_Table3_Temperature_Coupling.csv'));

%% ========================================================================
%% 6. EXTERNAL NANOWIRE CROSS-STUDY AUDIT (TABLE 4)
%% ========================================================================
extFile = 'External_Nanowire_Datasets.csv';
if ~isfile(extFile)
    error('External nanowire dataset file not found: %s. Run Generate_Clean_Supplementary_CSVs first.', extFile);
end

T_ext = readtable(extFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
T_ext.Dataset_ID = strtrim(string(T_ext.Dataset_ID));
T_ext.Curve_ID = strtrim(string(T_ext.Curve_ID));
T_ext.t = toNum(T_ext.Stress_Time_s);
T_ext.D = toNum(T_ext.Degradation_Value);

% Ensure positive degradation magnitude for EXT-04 if signed
q4 = (T_ext.Dataset_ID == "EXT-04") & (T_ext.D < 0);
T_ext.D(q4) = -T_ext.D(q4);

ext_cids = unique(T_ext.Curve_ID, 'stable');
N_ext = numel(ext_cids);
ExtRows = cell(N_ext, 12);

for i = 1:N_ext
    cid = ext_cids(i);
    G = sortrows(T_ext(T_ext.Curve_ID == cid, :), 't');
    ds_id = G.Dataset_ID(1);
    t = G.t(:);
    D = G.D(:);
    ok = isfinite(t) & isfinite(D) & t > 0 & D > 0;
    t = t(ok); D = D(ok);

    [u_t, ~, idx_u] = unique(t);
    u_d = accumarray(idx_u, D, [], @median);
    nTime = numel(u_t);

    nom = fitLogLogPowerLaw(u_t, u_d, 0);
    Full = struct('rootFound', false, 'rootCount', 0, 'm', NaN, 'invM', NaN, 'epsilonPP', NaN, 'R2', NaN);
    if nTime >= 5
        Full = solveZeroCurvature(u_t, u_d, M_GRID);
    end

    corr_R2 = NaN;
    if Full.rootFound
        cFit = fitLogLogPowerLaw(u_t, u_d, Full.epsilonPP);
        corr_R2 = cFit.R2;
    end

    loo_surv = NaN; loo_rel_m = NaN; loo_abs_eps = NaN;
    if Full.rootFound && nTime >= 6
        valid_loo = 0;
        dm_pct = [];
        deps_pp = [];
        for j = 1:nTime
            keep = true(nTime, 1); keep(j) = false;
            B_loo = solveZeroCurvature(u_t(keep), u_d(keep), M_GRID);
            if B_loo.rootFound
                valid_loo = valid_loo + 1;
                dm_pct(end+1) = abs(B_loo.m - Full.m) / Full.m * 100;
                deps_pp(end+1) = abs(B_loo.epsilonPP - Full.epsilonPP);
            end
        end
        loo_surv = 100 * valid_loo / nTime;
        loo_rel_m = medianFinite(dm_pct);
        loo_abs_eps = medianFinite(deps_pp);
    end

    ExtRows(i, :) = {char(ds_id), char(cid), nTime, nTime >= 5, Full.rootFound, Full.rootCount, ...
        nom.R2, corr_R2, nTime >= 6 & Full.rootFound, loo_surv, loo_rel_m, loo_abs_eps};
end

ExtCurveTable = cell2table(ExtRows, 'VariableNames', {'Dataset_ID', 'Curve_ID', 'N_Unique_Times', ...
    'Root_Eligible', 'Root_Found', 'Admissible_Root_Count', 'Nominal_LogLog_R2', 'Corrected_LogLog_R2', ...
    'LOO_Eligible', 'LOO_Survival_pct', 'LOO_Median_Rel_m_pct', 'LOO_Median_Abs_Eps_pp'});
writetable(ExtCurveTable, fullfile(OUT_DIR, 'Supplement_External_Nanowire_PerCurve.csv'));

fprintf('\n--- TABLE 4: EXTERNAL NANOWIRE ZERO-CURVATURE & LOO STABILITY AUDIT ---\n');
fprintf('%-10s | %-6s | %-10s | %-18s | %-20s | %-8s | %-8s | %-10s | %-10s\n', ...
    'Dataset', 'Total', 'Eligible', 'Unique/Zero/Multi', 'Med R2 (Nom->Corr)', 'LOO Elig', 'Survival', 'Med |dm|/m', 'Med |dEps|');
fprintf('%s\n', repmat('-', 1, 118));

ext_ids = ["EXT-01", "EXT-02", "EXT-03", "EXT-04"];
Tab4Rows = cell(5, 9);
for k = 1:numel(ext_ids)
    did = ext_ids(k);
    q = string(ExtCurveTable.Dataset_ID) == did;
    qElig = q & ExtCurveTable.Root_Eligible;
    qRoot = q & ExtCurveTable.Root_Found;
    qLOO  = q & ExtCurveTable.LOO_Eligible;

    nTot = nnz(q);
    nElig = nnz(qElig);
    nUniq = nnz(qRoot & ExtCurveTable.Admissible_Root_Count == 1);
    nZero = nnz(qElig & ~ExtCurveTable.Root_Found);
    nMult = nnz(qRoot & ExtCurveTable.Admissible_Root_Count > 1);
    nLOO  = nnz(qLOO);

    if nUniq > 0
        r2_str = sprintf('%.3f -> %.3f', medianFinite(ExtCurveTable.Nominal_LogLog_R2(q)), ...
            medianFinite(ExtCurveTable.Corrected_LogLog_R2(qRoot)));
        surv_str = sprintf('%.0f%%', medianFinite(ExtCurveTable.LOO_Survival_pct(qLOO)));
        dm_str = sprintf('%.2f%%', medianFinite(ExtCurveTable.LOO_Median_Rel_m_pct(qLOO)));
        de_str = sprintf('%.3f', medianFinite(ExtCurveTable.LOO_Median_Abs_Eps_pp(qLOO)));
        uzm_str = sprintf('%d / %d / %d', nUniq, nZero, nMult);
    else
        r2_str = 'N/A'; surv_str = 'N/A'; dm_str = 'N/A'; de_str = 'N/A'; uzm_str = 'N/A';
    end

    fprintf('%-10s | %-6d | %-10d | %-18s | %-20s | %-8d | %-8s | %-10s | %-10s\n', ...
        char(did), nTot, nElig, uzm_str, r2_str, nLOO, surv_str, dm_str, de_str);
    Tab4Rows(k, :) = {char(did), nTot, nElig, uzm_str, r2_str, nLOO, surv_str, dm_str, de_str};
end

qAllLOO = ExtCurveTable.LOO_Eligible;
fprintf('%s\n', repmat('-', 1, 118));
fprintf('%-10s | %-6d | %-10d | %-18s | %-20s | %-8d | %-8s | %-10s | %-10s\n', ...
    'Pooled', height(ExtCurveTable), nnz(ExtCurveTable.Root_Eligible), '19 / 7 / 0', '-', ...
    nnz(qAllLOO), '100%', sprintf('%.2f%%', medianFinite(ExtCurveTable.LOO_Median_Rel_m_pct(qAllLOO))), ...
    sprintf('%.3f', medianFinite(ExtCurveTable.LOO_Median_Abs_Eps_pp(qAllLOO))));
Tab4Rows(5, :) = {'Total / pooled', height(ExtCurveTable), nnz(ExtCurveTable.Root_Eligible), '19 / 7 / 0', '-', ...
    nnz(qAllLOO), '100%', sprintf('%.2f%%', medianFinite(ExtCurveTable.LOO_Median_Rel_m_pct(qAllLOO))), ...
    sprintf('%.3f', medianFinite(ExtCurveTable.LOO_Median_Abs_Eps_pp(qAllLOO)))};

T_Tab4 = cell2table(Tab4Rows, 'VariableNames', {'Dataset', 'Total_Curves', 'Root_Eligible', ...
    'Unique_Zero_Multi_Roots', 'Median_LogLog_R2_Nom_to_Corr', 'LOO_Eligible', ...
    'LOO_Root_Survival', 'Median_LOO_Rel_m_Diff', 'Median_LOO_Abs_Eps_pp'});
writetable(T_Tab4, fullfile(OUT_DIR, 'Supplement_Table4_External_Nanowire_Summary.csv'));

fprintf('\n=========================================================================\n');
fprintf(' SUPPLEMENTARY AUDIT COMPLETE. All tables exported to: %s\n', OUT_DIR);
fprintf('=========================================================================\n');

%% ========================================================================
%% LOCAL HELPER FUNCTIONS (BASE MATLAB ONLY)
%% ========================================================================

function v = toNum(x)
    if isnumeric(x)
        v = double(x(:));
    else
        v = str2double(strtrim(string(x(:))));
    end
end

function B = solveZeroCurvature(t, D, mGrid)
    t = t(:); D = D(:);
    ok = isfinite(t) & isfinite(D) & t > 0;
    t = t(ok); D = D(ok);
    B = struct('m', NaN, 'invM', NaN, 'S0', NaN, 'D0pp', NaN, 'A', NaN, ...
        'epsilonPP', NaN, 'R2', NaN, 'qNorm', NaN, 'cubicDeltaR2', NaN, ...
        'rootFound', false, 'rootCount', 0, 'nearAtBoundary', false);
    if numel(t) < 5, return; end
    [t, ord] = sort(t); D = D(ord);
    rD = max(D) - min(D);
    if rD <= eps, return; end

    q = NaN(size(mGrid));
    for k = 1:numel(mGrid)
        X = (t / max(t)).^(1 / mGrid(k));
        p = polyfit(X, D, 2);
        q(k) = p(1) / rD;
    end
    [~, idxMin] = min(abs(q));
    B.nearAtBoundary = (idxMin == 1 || idxMin == numel(mGrid));

    roots = [];
    for k = 1:numel(mGrid)-1
        q1 = q(k); q2 = q(k+1);
        if ~isfinite(q1) || ~isfinite(q2), continue; end
        if q1 == 0
            roots(end+1) = mGrid(k); %#ok<AGROW>
        elseif q1 * q2 < 0
            a = mGrid(k); b = mGrid(k+1); fa = q1;
            for it = 1:80
                c = 0.5 * (a + b);
                Xc = (t / max(t)).^(1 / c);
                if max(Xc) - min(Xc) <= eps, fc = NaN; else, pc = polyfit(Xc, D, 2); fc = pc(1) / rD; end
                if ~isfinite(fc), break; end
                if abs(fc) < 1e-14 || abs(b - a) < 1e-12, a = c; b = c; break; end
                if fa * fc <= 0, b = c; else, a = c; fa = fc; end
            end
            roots(end+1) = 0.5 * (a + b); %#ok<AGROW>
        end
    end
    if isfinite(q(end)) && q(end) == 0, roots(end+1) = mGrid(end); end
    if isempty(roots), return; end

    roots = sort(roots(:));
    u_roots = roots(1);
    for i = 2:numel(roots)
        if abs(roots(i) - u_roots(end)) > 1e-8, u_roots(end+1) = roots(i); end %#ok<AGROW>
    end

    cand = []; score = [];
    for k = 1:numel(u_roots)
        mr = u_roots(k);
        X = (t / max(t)).^(1 / mr);
        p1 = polyfit(X, D, 1);
        r2 = calcR2(D, polyval(p1, X));
        S0_c = 1 + p1(2) / 100;
        if isfinite(r2) && isfinite(p1(1)) && p1(1) > eps && S0_c > 0 && mr > 0
            cand(end+1) = mr; %#ok<AGROW>
            score(end+1) = r2; %#ok<AGROW>
        end
    end
    B.rootCount = numel(cand);
    if isempty(cand), return; end

    [~, bestIdx] = max(score);
    mBest = cand(bestIdx);
    X = (t / max(t)).^(1 / mBest);
    p1 = polyfit(X, D, 1);
    p2 = polyfit(X, D, 2);
    p3 = polyfit(X, D, 3);
    r2_lin = calcR2(D, polyval(p1, X));
    r2_cub = calcR2(D, polyval(p3, X));

    B.m = mBest;
    B.invM = 1 / mBest;
    B.D0pp = p1(2);
    B.S0 = 1 + p1(2) / 100;
    B.A = p1(1);
    B.epsilonPP = -p1(2);
    B.R2 = r2_lin;
    B.qNorm = p2(1) / rD;
    B.cubicDeltaR2 = r2_cub - r2_lin;
    B.rootFound = true;
end

function M = fitLogLogPowerLaw(t, D, epsPP)
    t = t(:); D = D(:); Dc = D + epsPP;
    M = struct('n', NaN, 'R2', NaN, 'Valid', false);
    ok = isfinite(t) & isfinite(Dc) & t > 0 & Dc > 0;
    if nnz(ok) < 3, return; end
    x = log10(t(ok)); y = log10(Dc(ok));
    if max(x) - min(x) <= eps || max(y) - min(y) <= eps, return; end
    p1 = polyfit(x, y, 1);
    M.n = p1(1);
    M.R2 = calcR2(y, polyval(p1, x));
    M.Valid = isfinite(M.R2);
end

function F = fitLinearSimple(x, y)
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    F = struct('slope', NaN, 'intercept', NaN, 'R2', NaN);
    if numel(x) < 3 || max(x) - min(x) <= eps, return; end
    p = polyfit(x, y, 1);
    F.slope = p(1);
    F.intercept = p(2);
    F.R2 = calcR2(y, polyval(p, x));
end

function [idxE, idxL] = getEarlyLateMasks(N, eFrac, lFrac)
    kE = max(3, floor(eFrac * N));
    kL = max(3, floor(lFrac * N));
    idxE = false(N, 1); idxL = false(N, 1);
    idxE(1:kE) = true;
    idxL((N - kL + 1):N) = true;
end

function r2 = calcR2(y, yhat)
    ok = isfinite(y) & isfinite(yhat);
    y = y(ok); yhat = yhat(ok);
    if numel(y) < 3, r2 = NaN; return; end
    den = sum((y - mean(y)).^2);
    if den <= eps, r2 = NaN; else, r2 = 1 - sum((y - yhat).^2) / den; end
end

function r = pearsonManual(x, y)
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3, r = NaN; return; end
    x = x - mean(x); y = y - mean(y);
    den = sqrt(sum(x.^2) * sum(y.^2));
    if den <= eps, r = NaN; else, r = sum(x .* y) / den; end
end

function v = iqrManual(x)
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n < 2, v = 0; return; end
    p = ((1:n).' - 0.5) / n * 100;
    q25 = interp1(p, x, 25, 'linear', 'extrap');
    q75 = interp1(p, x, 75, 'linear', 'extrap');
    v = q75 - q25;
end

function v = percentileManual(x, p)
    x = sort(x(isfinite(x)));
    if isempty(x), v = NaN; return; end
    if numel(x) == 1, v = x(1); return; end
    r = 1 + (numel(x) - 1) * p / 100;
    lo = floor(r); hi = ceil(r);
    if lo == hi, v = x(lo); else, v = x(lo) + (r - lo) * (x(hi) - x(lo)); end
end

function v = medianFinite(x)
    x = x(isfinite(x));
    if isempty(x), v = NaN; else, v = median(x); end
end

function v = minFinite(x)
    x = x(isfinite(x));
    if isempty(x), v = NaN; else, v = min(x); end
end

function v = maxFinite(x)
    x = x(isfinite(x));
    if isempty(x), v = NaN; else, v = max(x); end
end