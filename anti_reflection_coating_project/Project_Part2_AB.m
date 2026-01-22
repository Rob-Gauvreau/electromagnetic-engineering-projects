%% PART 2 (A) AND (B)  Double-Layer
clear; clc; close all;

n_0 = 1;
n_1 = 1.4;
n_2 = 2.62;
n_3 = 3.5;
j=1i;
Lambda_C = 650;
Lambda_Start = 200;
Lambda_End   = 2200;
Lambda_Range = Lambda_Start:Lambda_End;
g01 = (n_0 - n_1)./(n_0 + n_1);
g12 = (n_1 - n_2)./(n_1 + n_2);
g2S = (n_2 - n_3)./(n_2 + n_3);
t01 = 2*n_0./(n_0 + n_1);
t12 = 2*n_1./(n_1 + n_2);
t2S = 2*n_2./(n_2 + n_3);
Q01 = (1/t01)*([1 g01; g01 1]);
Q12 = (1/t12)*([1 g12; g12 1]);
Q2S = (1/t2S)*([1 g2S; g2S 1]);
Delta = (pi/2)*(Lambda_C./Lambda_Range);
Pcol  = [exp(j*Delta); exp(-j*Delta)];
Reflectance = zeros(size(Lambda_Range));
Power = zeros(size(Lambda_Range));
for k = 1:length(Lambda_Range)
    L = Lambda_Range(k);
    Pm = [Pcol(1,k) 0; 0 Pcol(2,k)];
    T = Q01 * Pm * Q12 * Pm * Q2S;
    r = T(2,1)./ T(1,1);
    R = abs(r).^2;
    Reflectance(k) = R;
    IRRAD = (6.16e15)/((L^5)*(exp(2484/L)-1));
    Power(k) = (1 - R)*IRRAD;
end

%% REFLECTIVITY GRAPH
Lambda_Start2 = 400;
Lambda_End2   = 1400;
idx2 = (Lambda_Range>=Lambda_Start2 & Lambda_Range<=Lambda_End2);

figure;
plot(Lambda_Range(idx2), Reflectance(idx2)*100,'LineWidth',2); grid on; hold on;
xline(Lambda_C,'r','LineWidth',4,'Label','\lambda_C = 650 nm', 'LabelOrientation','horizontal','FontSize',22,'LabelHorizontalAlignment','right');
title('Reflectivity vs Wavelength (400 nm - 1400 nm)','FontSize',22);
xlabel('Wavelength (nm)','FontSize',22);
ylabel('Reflectivity (%)','FontSize',22);
xlim([400,1600]);

%% POWER 200 - 2200
figure;
plot(Lambda_Range, Power,'LineWidth',2); grid on;
title('Transmitted Power vs Wavelength (200 nm - 2200 nm)','FontSize',22);
xlabel('Wavelength (nm)','FontSize',22);
ylabel('Power (W/m^2)','FontSize',22);
xlim([Lambda_Start,Lambda_End]);

%% POWER 400 - 1400
figure;
plot(Lambda_Range(idx2), Power(idx2),'LineWidth',2); grid on;
title('Transmitted Power vs Wavelength (400 nm - 1400 nm)','FontSize',22);
xlabel('Wavelength (nm)','FontSize',22);
ylabel('Power (W/m^2)','FontSize',22);
xlim([Lambda_Start2,Lambda_End2]);

%% TOTAL POWER PRINT
fprintf('Total Transmitted Power (200–2200 nm) = %.6f W/m^2\n', sum(Power));
fprintf('Total Transmitted Power (400–1400 nm) = %.6f W/m^2\n', sum(Power(idx2)));


