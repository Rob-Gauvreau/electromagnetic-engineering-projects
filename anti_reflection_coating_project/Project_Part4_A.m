%% Part 3 (A) AND (B) — Triple Layer
clc; clear; close all;
n_0 = 1.0;     
n_1 = 1.4;    
n_3 = 3.15;   
n_4 = 3.5;     
Lambda_C = 650;
j = 1i;

%% (A)
r01 = (n_0 - n_1) / (n_0 + n_1);
t01 = 2*n_0 / (n_0 + n_1);
Q01 = (1/t01) * [1 r01; r01 1];
r3S = (n_3 - n_4) / (n_3 + n_4);
t3S = 2*n_3 / (n_3 + n_4);
Q3S = (1/t3S) * [1 r3S; r3S 1];
Delta = pi/2;
P = [exp(j*Delta) 0; 0 exp(-j*Delta)];
n_2_range = 0:0.01:4.50;
Store_Reflectance = zeros(1, numel(n_2_range));
for i = 1:numel(n_2_range)
    n_2 = n_2_range(i);
    r12 = (n_1 - n_2) / (n_1 + n_2);
    t12 = 2*n_1 / (n_1 + n_2);
    Q12 = (1/t12) * [1 r12; r12 1];
    r23 = (n_2 - n_3) / (n_2 + n_3);
    t23 = 2*n_2 / (n_2 + n_3);
    Q23 = (1/t23) * [1 r23; r23 1];
    T = Q01 * P * Q12 * P * Q23 * P * Q3S;
    Gamma = T(2,1) / T(1,1);
    Store_Reflectance(i) = abs(Gamma)^2;
end
[~, Min_Index] = min(Store_Reflectance);
min_n_2 = n_2_range(Min_Index);
figure;
plot(n_2_range, Store_Reflectance * 100, 'LineWidth', 2);
grid on;
title('Reflectivity vs n_2 at \lambda_C = 650 nm (Triple Layer)', 'FontSize', 22);
xlabel('n_2', 'FontSize', 22);
ylabel('Reflectivity (%)', 'FontSize', 22);

fprintf('n_1 = %.2f, n_3 = %.2f, n_4 = %.2f\n', n_1, n_3, n_4);
fprintf('Minimum Reflectivity at \\lambda_C: n_2 = %.2f\n\n', min_n_2);

%% (B)
num_N_2 = 300;
Store_n_2 = linspace(1.4, 3.0, num_N_2); 
Store_Total_Power_400_1400 = computeTotalPower( ...
    400, 1400, Q01, Q3S, n_0, n_1, n_3, n_4, Lambda_C, num_N_2, Store_n_2);
Store_Total_Power_200_2200 = computeTotalPower( ...
    200, 2200, Q01, Q3S, n_0, n_1, n_3, n_4, Lambda_C, num_N_2, Store_n_2);
figure;
plot(Store_n_2, Store_Total_Power_400_1400, 'r', 'LineWidth', 1.8); hold on;
plot(Store_n_2, Store_Total_Power_200_2200, 'b', 'LineWidth', 1.8);
grid on;
title('Total Transmitted Power vs n_2 (Triple Layer)', 'FontSize', 22);
xlabel('n_2', 'FontSize', 22);
ylabel('Total Power (W/m^2)', 'FontSize', 22);
legend('Wavelength 400–1400 nm', 'Wavelength 200–2200 nm', 'Location', 'best');

[max_Power_1, idx1] = max(Store_Total_Power_400_1400);
fprintf('Max Power (400–1400 nm): %.2f W/m^2 at n_2 = %.2f\n', max_Power_1, Store_n_2(idx1));

[max_Power_2, idx2] = max(Store_Total_Power_200_2200);
fprintf('Max Power (200–2200 nm): %.2f W/m^2 at n_2 = %.2f\n', max_Power_2, Store_n_2(idx2));

%%  Helper
function Store_Total_Power = computeTotalPower(Lambda_Start, Lambda_End, Q01, Q3S, n_0, n_1, n_3, n_4, Lambda_C, numN2, Store_n_2)
    Lambda_Array = Lambda_Start:Lambda_End;
    Delta_Array  = (pi/2) * (Lambda_C ./ Lambda_Array);
    IRRAD_Array = 6.16e15 ./ ( Lambda_Array.^5 .* (exp(2484 ./ Lambda_Array) - 1) );
    Store_Total_Power = zeros(1, numN2);
    for idx = 1:numN2
        n_2 = Store_n_2(idx);
        r12 = (n_1 - n_2) / (n_1 + n_2);
        r23 = (n_2 - n_3) / (n_2 + n_3);
        t12 = 2*n_1 / (n_1 + n_2);
        t23 = 2*n_2 / (n_2 + n_3);
        Q12 = (1/t12) * [1 r12; r12 1];
        Q23 = (1/t23) * [1 r23; r23 1];
        Pcol = [exp(1j * Delta_Array); exp(-1j * Delta_Array)];
        Tcol = zeros(2, numel(Lambda_Array)); % columns hold field at input
        for k = 1:numel(Lambda_Array)
            P = [Pcol(1, k) 0; 0 Pcol(2, k)];
            Tcol(:, k) = Q01 * P * Q12 * P * Q23 * P * Q3S * [1; 0];
        end
        Tau   = 1 ./ Tcol(1, :);
        Trans = (abs(Tau).^2) * (n_4 / n_0);
        Store_Total_Power(idx) = sum(Trans .* IRRAD_Array);
    end
end