%% PART 2 (c)  – Optimize n1 and n2 (double layer)
clear; clc; close all;
n_0 = 1;
n_3 = 3.5;
Lambda_C = 650;
j = 1i;

IRRAD_func = @(L) (6.16e15)./((L.^5) .* (exp(2484./L) - 1));
Phase_func = @(Lambda_vec) (pi/2) * (Lambda_C ./ Lambda_vec);
n1_min0 = 1.0;  n1_max0 = 3.0;
n2_min0 = 1.0;  n2_max0 = 3.0;
step0   = 0.4;
Max_Iteration = 5;
min_step = 0.01;

%% OPTIMIZATION #1: 200–2200 nm
Lambda_Start_A = 200;
Lambda_End_A   = 2200;
Lambda_Range_A = Lambda_Start_A:Lambda_End_A;
n1_min = n1_min0; n1_max = n1_max0;
n2_min = n2_min0; n2_max = n2_max0;
step   = step0;
bestPower_A = -Inf; best_n1_A = NaN; best_n2_A = NaN; bestR_A = [];
for iter = 1:Max_Iteration
    n1_vals = n1_min:step:n1_max;
    n2_vals = n2_min:step:n2_max;
    localBestP = -Inf; localBest_n1 = NaN; localBest_n2 = NaN; localBestR = [];
    for n_1 = n1_vals
        for n_2 = n2_vals
            [R_vec, P_vec] = compute_reflect_trans_power( ...
                n_0, n_1, n_2, n_3, Lambda_C, Lambda_Range_A, IRRAD_func, Phase_func, j);
            Psum = sum(P_vec);
            if Psum > localBestP
                localBestP = Psum;
                localBest_n1 = n_1;
                localBest_n2 = n_2;
                localBestR = R_vec;
            end
        end
    end
    if localBestP > bestPower_A
        bestPower_A = localBestP;
        best_n1_A = localBest_n1;
        best_n2_A = localBest_n2;
        bestR_A   = localBestR;
    end
    n1_min = max(localBest_n1 - 2*step, n1_min0);
    n1_max = min(localBest_n1 + 2*step, n1_max0);
    n2_min = max(localBest_n2 - 2*step, n2_min0);
    n2_max = min(localBest_n2 + 2*step, n2_max0);
    step = max(step/2, min_step);
end

[bestRcurve_A, bestPcurve_A] = compute_reflect_trans_power(n_0, best_n1_A, best_n2_A, n_3, Lambda_C, Lambda_Range_A, IRRAD_func, Phase_func, j);
bestPtot_A = sum(bestPcurve_A);

%% OPTIMIZATION #2: 400–1400 nm
Lambda_Start_B = 400;
Lambda_End_B   = 1400;
Lambda_Range_B = Lambda_Start_B:Lambda_End_B;
n1_min = n1_min0; n1_max = n1_max0;
n2_min = n2_min0; n2_max = n2_max0;
step   = step0;
bestPower_B = -Inf; best_n1_B = NaN; best_n2_B = NaN; bestR_B = [];
for iter = 1:Max_Iteration
    n1_vals = n1_min:step:n1_max;
    n2_vals = n2_min:step:n2_max;
    localBestP = -Inf; localBest_n1 = NaN; localBest_n2 = NaN; localBestR = [];
    for n_1 = n1_vals
        for n_2 = n2_vals
            [R_vec, P_vec] = compute_reflect_trans_power( ...
                n_0, n_1, n_2, n_3, Lambda_C, Lambda_Range_B, IRRAD_func, Phase_func, j);
            Psum = sum(P_vec);
            if Psum > localBestP
                localBestP = Psum;
                localBest_n1 = n_1;
                localBest_n2 = n_2;
                localBestR = R_vec;
            end
        end
    end
    if localBestP > bestPower_B
        bestPower_B = localBestP;
        best_n1_B = localBest_n1;
        best_n2_B = localBest_n2;
        bestR_B   = localBestR;
    end
    n1_min = max(localBest_n1 - 2*step, n1_min0);
    n1_max = min(localBest_n1 + 2*step, n1_max0);
    n2_min = max(localBest_n2 - 2*step, n2_min0);
    n2_max = min(localBest_n2 + 2*step, n2_max0);
    step = max(step/2, min_step);
end
[bestRcurve_B, bestPcurve_B] = compute_reflect_trans_power( ...
    n_0, best_n1_B, best_n2_B, n_3, Lambda_C, Lambda_Range_B, IRRAD_func, Phase_func, j);
bestPtot_B = sum(bestPcurve_B);

figure;
plot(Lambda_Range_A, bestRcurve_A*100, 'LineWidth', 2); grid on;
title('Optimized Reflectivity vs Wavelength (200 nm - 2200 nm)','FontSize',22);
xlabel('Wavelength (nm)','FontSize',22);
ylabel('Reflectivity (%)','FontSize',22);
xlim([Lambda_Start_A, Lambda_End_A]);

figure;
plot(Lambda_Range_B, bestRcurve_B*100, 'LineWidth', 2); grid on;
title('Optimized Reflectivity vs Wavelength (400 nm - 1400 nm)','FontSize',22);
xlabel('Wavelength (nm)','FontSize',22);
ylabel('Reflectivity (%)','FontSize',22);
xlim([Lambda_Start_B, Lambda_End_B]);

fprintf('=== Optimization over 200–2200 nm ===\n');
fprintf('Optimal n_1 = %.4f\n', best_n1_A);
fprintf('Optimal n_2 = %.4f\n', best_n2_A);
fprintf('Total Transmitted Power (200–2200 nm) = %.6f W/m^2\n\n', bestPtot_A);
fprintf('=== Optimization over 400–1400 nm ===\n');
fprintf('Optimal n_1 = %.4f\n', best_n1_B);
fprintf('Optimal n_2 = %.4f\n', best_n2_B);
fprintf('Total Transmitted Power (400–1400 nm) = %.6f W/m^2\n', bestPtot_B);

%% Helper Function
function [Reflectance_vec, Power_vec] = compute_reflect_trans_power( ...
    n_0, n_1, n_2, n_3, Lambda_C, Lambda_Range, IRRAD_func, Phase_func, j)
    g01 = (n_0 - n_1) / (n_0 + n_1);
    g12 = (n_1 - n_2) / (n_1 + n_2);
    g2S = (n_2 - n_3) / (n_2 + n_3);
    t01 = 2*n_0 / (n_0 + n_1);
    t12 = 2*n_1 / (n_1 + n_2);
    t2S = 2*n_2 / (n_2 + n_3);
    Q01 = (1/t01) * [1 g01; g01 1];
    Q12 = (1/t12) * [1 g12; g12 1];
    Q2S = (1/t2S) * [1 g2S; g2S 1];
    Delta = Phase_func(Lambda_Range);
    Pcol  = [exp(j*Delta); exp(-j*Delta)];
    Reflectance_vec = zeros(size(Lambda_Range));
    Power_vec       = zeros(size(Lambda_Range));
    for k = 1:length(Lambda_Range)
        L = Lambda_Range(k);
        Pm = [Pcol(1,k) 0; 0 Pcol(2,k)];
        T = Q01 * Pm * Q12 * Pm * Q2S;
        Gamma = T(2,1) / T(1,1);
        R = abs(Gamma)^2;
        Reflectance_vec(k) = R;
        Tau = 1 / T(1,1);
        Trans = abs(Tau)^2 * (n_3 / n_0);
        IRRAD = IRRAD_func(L);
        Power_vec(k) = Trans * IRRAD;
    end
end