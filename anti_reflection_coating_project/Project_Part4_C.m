%% Part 3 (c) – Triple Layer Optimization

%% Fixed indices and center wavelength
n0 = 1.0;    % air
n4 = 3.5;    % solar cell substrate
Lambda_C = 650; % nm
j = 1i;
IRRAD_fun = @(L) 6.16e15 ./ ( (L.^5) .* (exp(2484./L) - 1) );
phase_fun = @(Lvec) (pi/2) * (Lambda_C  ./ Lvec);   % Delta(λ) for quarter-wave scaling
n1_min0 = 1.0;  n1_max0 = 2.0;
n2_min0 = 1.5;  n2_max0 = 3.0;
n3_min0 = 2.0;  n3_max0 = 3.6;

step0   = 0.5;
min_step = 0.01;
MaxIteration = 5;

%% OPTIMIZATION #1: 200–2200 nm
LminA = 200;  LmaxA = 2200;
LvecA = LminA:LmaxA;
DeltaA = phase_fun(LvecA);                 % 1 x NA
pA = exp(j*DeltaA);                         % 1 x NA
invpA = exp(-j*DeltaA);
IRRAA = IRRAD_fun(LvecA);                  % 1 x NA

[best_n1_A, best_n2_A, best_n3_A, RbestA, PtotA] = optimize_triple(n0, n4, LvecA, pA, invpA, IRRAA, n1_min0, n1_max0, n2_min0, n2_max0, n3_min0, n3_max0, step0, min_step, MaxIteration);
figure;
plot(LvecA, RbestA*100, 'b', 'LineWidth', 2);
grid on;
title('Optimized Reflectivity vs Wavelength (200 nm - 2200 nm)', 'FontSize', 22);
xlabel('Wavelength (nm)', 'FontSize', 22);
ylabel('Reflectivity (%)', 'FontSize', 22);
xlim([LminA, LmaxA]);
fprintf('Triple Layer Optimization over 200–2200 nm\n');
fprintf('Best n1 = %.4f\n', best_n1_A);
fprintf('Best n2 = %.4f\n', best_n2_A);
fprintf('Best n3 = %.4f\n', best_n3_A);
fprintf('Total Transmitted Power (200–2200 nm) = %.6f W/m^2\n\n', PtotA);

%% OPTIMIZATION #2: 400–1400 nm 
LminB = 400;  LmaxB = 1400;
LvecB = LminB:LmaxB;
DeltaB = phase_fun(LvecB);
pB = exp(j*DeltaB);
invpB = exp(-j*DeltaB);
IRRAB = IRRAD_fun(LvecB);
[best_n1_B, best_n2_B, best_n3_B, RbestB, PtotB] = optimize_triple(n0, n4, LvecB, pB, invpB, IRRAB, n1_min0, n1_max0, n2_min0, n2_max0, n3_min0, n3_max0, step0, min_step, MaxIteration);
figure;
plot(LvecB, RbestB*100, 'b', 'LineWidth', 2);
grid on;
title('Optimized Reflectivity vs Wavelength (400 nm - 1400 nm)', 'FontSize', 22);
xlabel('Wavelength (nm)', 'FontSize', 22);
ylabel('Reflectivity (%)', 'FontSize', 22);
xlim([LminB, LmaxB]);
fprintf('Triple Layer Optimization over 400–1400 nm\n');
fprintf('Best n1 = %.4f\n', best_n1_B);
fprintf('Best n2 = %.4f\n', best_n2_B);
fprintf('Best n3 = %.4f\n', best_n3_B);
fprintf('Total Transmitted Power (400–1400 nm) = %.6f W/m^2\n', PtotB);

%% Helper
function [best_n1, best_n2, best_n3, R_best_curve, P_best_total] = optimize_triple(n0, n4, Lvec, p, invp, IRRAD_vec, n1_min0, n1_max0, n2_min0, n2_max0, n3_min0, n3_max0, step0, min_step, MaxIteration)
    step   = step0;
    n1_min = n1_min0;  n1_max = n1_max0;
    n2_min = n2_min0;  n2_max = n2_max0;
    n3_min = n3_min0;  n3_max = n3_max0;
    bestPower = -Inf;
    best_n1 = NaN; best_n2 = NaN; best_n3 = NaN;
    R_best_curve = [];
    for it = 1:MaxIteration
        n1_vals = n1_min:step:n1_max;
        n2_vals = n2_min:step:n2_max;
        n3_vals = n3_min:step:n3_max;
        localBest = -Inf;
        local_n1 = NaN; local_n2 = NaN; local_n3 = NaN;
        local_R_curve = [];
        for n1 = n1_vals
            for n2 = n2_vals
                for n3 = n3_vals
                    [R_curve, P_curve] = compute_triple_RP_vec( ...
                        n0, n1, n2, n3, n4, p, invp, IRRAD_vec);
                    Psum = sum(P_curve);
                    if Psum > localBest
                        localBest = Psum;
                        local_n1 = n1; local_n2 = n2; local_n3 = n3;
                        local_R_curve = R_curve;
                    end
                end
            end
        end
        if localBest > bestPower
            bestPower = localBest;
            best_n1 = local_n1; best_n2 = local_n2; best_n3 = local_n3;
            R_best_curve = local_R_curve;
        end
        n1_min = max(local_n1 - 2*step, n1_min0);
        n1_max = min(local_n1 + 2*step, n1_max0);
        n2_min = max(local_n2 - 2*step, n2_min0);
        n2_max = min(local_n2 + 2*step, n2_max0);
        n3_min = max(local_n3 - 2*step, n3_min0);
        n3_max = min(local_n3 + 2*step, n3_max0);
        step   = max(step/2, min_step);
    end
    [R_best_curve, P_curve] = compute_triple_RP_vec( ...
        n0, best_n1, best_n2, best_n3, n4, p, invp, IRRAD_vec);
    P_best_total = sum(P_curve);
end
%% Helper
function [R_vec, P_vec] = compute_triple_RP_vec(n0, n1, n2, n3, n4, p, invp, IRRAD)
    g01 = (n0 - n1) / (n0 + n1);  t01 = 2*n0 / (n0 + n1);
    g12 = (n1 - n2) / (n1 + n2);  t12 = 2*n1 / (n1 + n2);
    g23 = (n2 - n3) / (n2 + n3);  t23 = 2*n2 / (n2 + n3);
    g3S = (n3 - n4) / (n3 + n4);  t3S = 2*n3 / (n3 + n4);
    a1 = 1/t01; b1 = g01/t01;   
    a2 = 1/t12; b2 = g12/t12;   
    a3 = 1/t23; b3 = g23/t23;   
    a4 = 1/t3S; b4 = g3S/t3S;   
    m11 = a1; m12 = b1; m21 = b1; m22 = a1;
    c11 = m11 .* p;      c12 = m12 .* invp;
    c21 = m21 .* p;      c22 = m22 .* invp;

    d11 = c11.*a2 + c12.*b2;
    d12 = c11.*b2 + c12.*a2;
    d21 = c21.*a2 + c22.*b2;
    d22 = c21.*b2 + c22.*a2;
    e11 = d11 .* p;      e12 = d12 .* invp;
    e21 = d21 .* p;      e22 = d22 .* invp;
    f11 = e11.*a3 + e12.*b3;
    f12 = e11.*b3 + e12.*a3;
    f21 = e21.*a3 + e22.*b3;
    f22 = e21.*b3 + e22.*a3;
    g11 = f11 .* p;      g12 = f12 .* invp;
    g21 = f21 .* p;      g22 = f22 .* invp;
    T11 = g11.*a4 + g12.*b4;
    T12 = g11.*b4 + g12.*a4;
    T21 = g21.*a4 + g22.*b4;
    Gamma = T21 ./ T11;
    R_vec = abs(Gamma).^2;
    Tau = 1 ./ T11;
    Trans = abs(Tau).^2 * (n4 / n0);
    P_vec = Trans .* IRRAD;
end