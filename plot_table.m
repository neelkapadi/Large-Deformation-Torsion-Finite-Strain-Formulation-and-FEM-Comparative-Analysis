% =========================================================================
% MATLAB Script to Process and Plot Shaft Torsion Data
% Ensure 'shaft_torsion_data.csv' is in the same folder as this script.
% =========================================================================

clear; clc; close all;

% 1. Load Data
% 'PreserveVariableNames' ensures MATLAB doesn't alter your C++ column headers
disp('Loading data from shaft_torsion_data.csv...');
data = readtable('shaft_torsion_data.csv', 'PreserveVariableNames', true);

% Average the axisymmetric data across all Theta angles for each Radius (R)
avg_data = groupsummary(data, 'R', 'mean', ...
    {'F_Error_Frobenius', 'Tau_Small', 'Tau_Finite', 'Tau_Classical', ...
     'Sigma_ZZ_Small', 'Sigma_ZZ_Finite'});

radius = avg_data.R;

% =========================================================================
% PLOT 1: Error between Numerical and Analytical Result
% =========================================================================
figure('Name', 'F Matrix Error', 'Color', 'w');
plot(radius, avg_data.mean_F_Error_Frobenius, 'r-', 'LineWidth', 2);
title('Error between Numerical and Analytical Deformation Gradient');
xlabel('Radius R (m)');
ylabel('Frobenius Norm of Error');
grid on;

% Create Table for Report
disp('--- TABLE 1: Error vs Radius ---');
Table_Error = table(radius, avg_data.mean_F_Error_Frobenius, ...
    'VariableNames', {'Radius_m', 'F_Error_Frobenius'});
disp(head(Table_Error, 5)); % Display first 5 rows

% =========================================================================
% PLOT 2: Torsional Shear Stress (sigma_z_theta) vs Radius
% =========================================================================
figure('Name', 'Shear Stress vs Radius', 'Color', 'w');
hold on;
plot(radius, avg_data.mean_Tau_Small / 1e6, '--', 'LineWidth', 1.5, 'DisplayName', 'Small Strain');
plot(radius, avg_data.mean_Tau_Finite / 1e6, '-', 'LineWidth', 2, 'DisplayName', 'Finite Strain');
plot(radius, avg_data.mean_Tau_Classical / 1e6, ':', 'LineWidth', 1.5, 'DisplayName', 'Classical');
hold off;
title('Torsional Shear Stress (\tau_{z\theta}) vs. Radius');
xlabel('Radius R (m)');
ylabel('Shear Stress (MPa)');
legend('Location', 'northwest');
grid on;

% Create Table for Report
disp('--- TABLE 2: Shear Stress (MPa) vs Radius ---');
Table_Shear = table(radius, avg_data.mean_Tau_Small/1e6, avg_data.mean_Tau_Finite/1e6, avg_data.mean_Tau_Classical/1e6, ...
    'VariableNames', {'Radius_m', 'Small_Strain_MPa', 'Finite_Strain_MPa', 'Classical_MPa'});
disp(head(Table_Shear, 5));

% =========================================================================
% PLOT 3: Axial Normal Stress (sigma_zz) vs Radius
% =========================================================================
figure('Name', 'Axial Stress vs Radius', 'Color', 'w');
hold on;
plot(radius, avg_data.mean_Sigma_ZZ_Small / 1e6, '--', 'LineWidth', 1.5, 'DisplayName', 'Small Strain');
plot(radius, avg_data.mean_Sigma_ZZ_Finite / 1e6, '-', 'LineWidth', 2, 'DisplayName', 'Finite Strain');
hold off;
title('Axial Normal Stress (\sigma_{zz}) vs. Radius (Poynting Effect)');
xlabel('Radius R (m)');
ylabel('Axial Stress (MPa)');
legend('Location', 'northwest');
grid on;

% Create Table for Report
disp('--- TABLE 3: Axial Stress (MPa) vs Radius ---');
Table_Axial = table(radius, avg_data.mean_Sigma_ZZ_Small/1e6, avg_data.mean_Sigma_ZZ_Finite/1e6, ...
    'VariableNames', {'Radius_m', 'Small_Strain_MPa', 'Finite_Strain_MPa'});
disp(head(Table_Axial, 5));

% =========================================================================
% PLOT 4: Torque versus Twist (alpha)
% =========================================================================
% NOTE: Replace these array values with the manual runs from your C++ code
% over different input twist rates (alpha).
alphas        = [0.0, 0.5, 1.0, 1.5, 2.0];          % Twist rates (rad/m)
torque_small  = [0, 377123, 754247, 1.13137e+06, 1.50849e+06];  % Torque values (N.m)
torque_finite = [0, 377398, 756442, 1.13878e+06, 1.52606e+06];  
torque_class  = [0, 377595, 755191, 1.13279e+06, 1.52606e+06];  

figure('Name', 'Torque vs Twist', 'Color', 'w');
hold on;
plot(alphas, torque_small, 'o--', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Small Strain');
plot(alphas, torque_finite, 's-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Finite Strain');
plot(alphas, torque_class, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Classical');
hold off;
title('Torque vs. Rate of Twist');
xlabel('Rate of Twist, \alpha (rad/m)');
ylabel('Torque (N\cdotm)');
legend('Location', 'northwest');
grid on;

% Create Table for Report
disp('--- TABLE 4: Torque vs Twist ---');
Table_Torque = table(alphas', torque_small', torque_finite', torque_class', ...
    'VariableNames', {'Alpha_rad_m', 'Torque_Small_Nm', 'Torque_Finite_Nm', 'Torque_Class_Nm'});
disp(Table_Torque);

disp('All plots generated successfully!');