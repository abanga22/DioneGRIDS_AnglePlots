%% Housekeeping
clearvars;
clc;

%% Define w and d
w = 1.8;
d = 0.45;
C = w/(2*d);

%% Get all possible nominal current combinations
I_nominal = [90 92 94 96 98 100 102 104 106 108 110].*10^(-9);
combinations = combvec(I_nominal, I_nominal, I_nominal, I_nominal); % Generate current combinations

%% Get alpha and beta for each combination of four nominal currents.
alpha_nom = zeros(1,length(combinations(1,:))); % initialize a vector to store alpha
alpha_currents_nom = zeros(4,length(combinations(1,:))); % vector to hold current combinations that yield particular alphas
beta_nom = zeros(1,length(combinations(1,:))); % initialize a vector to store beta
beta_currents_nom = zeros(4,length(combinations(1,:))); % vector to hold current combinations that yield particular betas
for i = 1:length(combinations(1,:)) % loop through the current combinations
    currents = combinations(:,i); % store the currents for a particular combination
    I_a = currents(1);
    I_b = currents(2);
    I_c = currents(3);
    I_d = currents(4);
    alpha_nom(i) = atand(C*((I_c + I_d - (I_a + I_b))/(I_c + I_d + I_a + I_b))); % calculate/store alpha
    alpha_currents_nom(:,i) = [I_a; I_b; I_c; I_d]; % store the current combination
    beta_nom(i) = atand(C*((I_b + I_c - (I_a + I_d))/(I_b + I_c + I_a + I_d))); % calculate/store beta
    beta_currents_nom(:,i) = alpha_currents_nom(:,i); % store the current combination
end

%% Get medians from measured cal data

% Channel A
I_A_filename = "GRIDS_DIONE_CALIBRATION_12132022_1_A_90-110nA_L.xlsx";
I_A_sheet = "IDM Data (arrays)";
I_A_matrix = readmatrix(I_A_filename, "Sheet", I_A_sheet);
I_A = I_A_matrix(:,6);

I_A_90med = median(I_A(1:1288));
I_A_92med = median(I_A(1509:2632));
I_A_94med = median(I_A(2874:3860));
I_A_96med = median(I_A(4077:5036));
I_A_98med = median(I_A(5221:6685));
I_A_100med = median(I_A(6984:8039));
I_A_102med = median(I_A(8377:9188));
I_A_104med = median(I_A(9369:10200));
I_A_106med = median(I_A(10410:11344));
I_A_108med = median(I_A(11513:12324));
I_A_110med = median(I_A(12545:13320));

I_A_med = [I_A_90med I_A_92med I_A_94med I_A_96med I_A_98med I_A_100med I_A_102med I_A_104med I_A_106med I_A_108med I_A_110med];

% Channel B
I_B_filename = "GRIDS_DIONE_CALIBRATION_12142022_2_B_90-110nA_L.xlsx";
I_B_sheet = "IDM Data (arrays)";
I_B_matrix = readmatrix(I_B_filename, "Sheet", I_B_sheet);
I_B = I_B_matrix(:,7);

I_B_90med = median(I_B(1:1076));
I_B_92med = median(I_B(1206:2124));
I_B_94med = median(I_B(2309:3104));
I_B_96med = median(I_B(3309:4056));
I_B_98med = median(I_B(4317:5036));
I_B_100med = median(I_B(5297:5964));
I_B_102med = median(I_B(6205:7040));
I_B_104med = median(I_B(7197:7888));
I_B_106med = median(I_B(8065:8804));
I_B_108med = median(I_B(9007:9653));
I_B_110med = median(I_B(9853:10480));

I_B_med = [I_B_90med I_B_92med I_B_94med I_B_96med I_B_98med I_B_100med I_B_102med I_B_104med I_B_106med I_B_108med I_B_110med];

% Channel C
I_C_filename = "GRIDS_DIONE_CALIBRATION_12082022_3_C_90nA-110nA_L.xlsx";
I_C_sheet = "IDM Data (arrays)";
I_C_matrix = readmatrix(I_C_filename, "Sheet", I_C_sheet);
I_C = I_C_matrix(:,8);

I_C_90med = median(I_C(65:1208));
I_C_92med = median(I_C(1229:2336));
I_C_94med = median(I_C(2345:3466));
I_C_96med = median(I_C(3489:4504));
I_C_98med = median(I_C(4527:5366));
I_C_100med = median(I_C(5397:6448));
I_C_102med = median(I_C(6469:7416));
I_C_104med = median(I_C(7425:8731));
I_C_106med = median(I_C(8749:9800));
I_C_108med = median(I_C(9817:10788));
I_C_110med = median(I_C(10789:11160));

I_C_med = [I_C_90med I_C_92med I_C_94med I_C_96med I_C_98med I_C_100med I_C_102med I_C_104med I_C_106med I_C_108med I_C_110med];

% Channel D
I_D_filename = "GRIDS_DIONE_CALIBRATION_12142022_4_D_90-110nA_L.xlsx";
I_D_sheet = "IDM Data (arrays)";
I_D_matrix = readmatrix(I_D_filename, "Sheet", I_D_sheet);
I_D = I_D_matrix(:,9);

I_D_90med = median(I_D(33:740));
I_D_92med = median(I_D(929:1564));
I_D_94med = median(I_D(1805:2696));
I_D_96med = median(I_D(2906:3548));
I_D_98med = median(I_D(3769:4484));
I_D_100med = median(I_D(4713:5468));
I_D_102med = median(I_D(5661:6300));
I_D_104med = median(I_D(6473:7140));
I_D_106med = median(I_D(7361:7996));
I_D_108med = median(I_D(8188:8772));
I_D_110med = median(I_D(8961:9600));

I_D_med = [I_D_90med I_D_92med I_D_94med I_D_96med I_D_98med I_D_100med I_D_102med I_D_104med I_D_106med I_D_108med I_D_110med];

%% Replace current combination values with appropriate measured currents
combinations_measured = combinations;

% Collector A
for i = 1:length(combinations(1,:))
    if combinations(1,i) - 90/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(1);
    elseif combinations(1,i) - 92/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(2);
    elseif combinations(1,i) - 94/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(3);
    elseif combinations(1,i) - 96/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(4);
    elseif combinations(1,i) - 98/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(5);
    elseif combinations(1,i) - 100/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(6);
    elseif combinations(1,i) - 102/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(7);
    elseif combinations(1,i) - 104/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(8);
    elseif combinations(1,i) - 106/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(9);
    elseif combinations(1,i) - 108/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(10);
    elseif combinations(1,i) - 110/10^9 <= 10^(-20)
        combinations_measured(1,i) = I_A_med(11);
    end
end

% Collector B
for i = 1:length(combinations(2,:))
    if combinations(2,i) - 90/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(1);
    elseif combinations(2,i) - 92/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(2);
    elseif combinations(2,i) - 94/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(3);
    elseif combinations(2,i) - 96/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(4);
    elseif combinations(2,i) - 98/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(5);
    elseif combinations(2,i) - 100/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(6);
    elseif combinations(2,i) - 102/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(7);
    elseif combinations(2,i) - 104/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(8);
    elseif combinations(2,i) - 106/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(9);
    elseif combinations(2,i) - 108/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(10);
    elseif combinations(2,i) - 110/10^9 <= 10^(-20)
        combinations_measured(2,i) = I_B_med(11);
    end
end

% Collector C
for i = 1:length(combinations(3,:))
    if combinations(3,i) - 90/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(1);
    elseif combinations(3,i) - 92/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(2);
    elseif combinations(3,i) - 94/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(3);
    elseif combinations(3,i) - 96/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(4);
    elseif combinations(3,i) - 98/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(5);
    elseif combinations(3,i) - 100/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(6);
    elseif combinations(3,i) - 102/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(7);
    elseif combinations(3,i) - 104/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(8);
    elseif combinations(3,i) - 106/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(9);
    elseif combinations(3,i) - 108/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(10);
    elseif combinations(3,i) - 110/10^9 <= 10^(-20)
        combinations_measured(3,i) = I_C_med(11);
    end
end

% Collector D
for i = 1:length(combinations(4,:))
    if combinations(4,i) - 90/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(1);
    elseif combinations(4,i) - 92/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(2);
    elseif combinations(4,i) - 94/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(3);
    elseif combinations(4,i) - 96/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(4);
    elseif combinations(4,i) - 98/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(5);
    elseif combinations(4,i) - 100/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(6);
    elseif combinations(4,i) - 102/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(7);
    elseif combinations(4,i) - 104/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(8);
    elseif combinations(4,i) - 106/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(9);
    elseif combinations(4,i) - 108/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(10);
    elseif combinations(4,i) - 110/10^9 <= 10^(-20)
        combinations_measured(4,i) = I_D_med(11);
    end
end

%% Get alpha and beta for each combination of four measured currents.
alpha_measured = zeros(1,length(combinations(1,:))); % initialize a vector to store alpha
alpha_currents_measured = zeros(4,length(combinations(1,:))); % vector to hold current combinations that yield particular alphas
beta_measured = zeros(1,length(combinations(1,:))); % initialize a vector to store beta
beta_currents_measured = zeros(4,length(combinations(1,:))); % vector to hold current combinations that yield particular betas
for i = 1:length(combinations_measured(1,:)) % loop through the current combinations
    currents = combinations_measured(:,i); % store the currents for a particular combination
    I_a = currents(1);
    I_b = currents(2);
    I_c = currents(3);
    I_d = currents(4);
    alpha_measured(i) = atand(C*((I_c + I_d - (I_a + I_b))/(I_c + I_d + I_a + I_b))); % calculate/store alpha
    alpha_currents_measured(:,i) = [I_a; I_b; I_c; I_d]; % store the current combination
    beta_measured(i) = atand(C*((I_b + I_c - (I_a + I_d))/(I_b + I_c + I_a + I_d))); % calculate/store beta
    beta_currents_measured(:,i) = alpha_currents_measured(:,i); % store the current combination
end

%% Plot raw data

figure(1);
plot(sort(alpha_nom), sort(alpha_nom), 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
scatter(sort(alpha_nom), sort(alpha_measured), 'LineWidth', 2);
title("Nominal vs. Measured α (raw data)");
xlabel("Nominal Angle (Degrees)");
ylabel("Measured Angle (Degrees)");
xlim([-5 5]);
ylim([-5 5]);
legend("Ideal Case", "Real Case");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

figure(2);
plot(sort(beta_nom), sort(beta_nom), 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
scatter(sort(beta_nom), sort(beta_measured));
title("Nominal vs. Measured β (raw data)");
xlabel("Nominal Angle (Degrees)");
ylabel("Measured Angle (Degrees)");
xlim([-5 5]);
ylim([-5 5]);
legend("Ideal Case", "Real Case");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

%% Plotting the data raw, we can see that for a given nominal angle, there 
%% are multiple measured angle possibilities. This is because there are 
%% multiple nominal current combinations that yield the same angle. To correct
%% for this, we take the median of the measured currents at a given nominal current.

alpha_measured_med = zeros(1,length(alpha_nom));
alpha_maxes = zeros(1,length(alpha_nom));
alpha_mins = zeros(1,length(alpha_nom));
for i = 1:length(alpha_nom)
    alpha_arr = alpha_measured(abs(alpha_nom - alpha_nom(i)) < 1e-10);
    alpha_med = median(alpha_arr);
    alpha_maxes(i) = max(alpha_arr);
    alpha_mins(i) = min(alpha_arr);
    alpha_measured_med(i) = alpha_med;
end

beta_measured_med = zeros(1,length(beta_nom));
beta_maxes = zeros(1,length(alpha_nom));
beta_mins = zeros(1,length(alpha_nom));
for i = 1:length(beta_nom)
    beta_arr = beta_measured(abs(beta_nom - beta_nom(i)) < 1e-10);
    beta_med = median(beta_arr);
    beta_maxes(i) = max(beta_arr);
    beta_mins(i) = min(beta_arr);
    beta_measured_med(i) = beta_med;
end

%% Plot median data

% Correct for alpha error if needed
% alpha_measured_med = alpha_measured_med + 0.0264;

figure(3);
plot(sort(alpha_nom), sort(alpha_nom), 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
errorbar(sort(alpha_nom), sort(alpha_measured_med), sort(alpha_measured_med) - sort(alpha_mins), sort(alpha_maxes) - sort(alpha_measured_med), 'o-');
title("Nominal vs. Measured α IDM Median");
xlabel("Nominal Angle (Degrees)");
ylabel("Measured Angle (Degrees)");
xlim([-5 5]);
ylim([-5 5]);
legend("Ideal Case", "Real Case");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

figure(4);
plot(sort(beta_nom), sort(beta_nom), 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
errorbar(sort(beta_nom), sort(beta_measured_med), sort(beta_measured_med) - sort(beta_mins), sort(beta_maxes) - sort(beta_measured_med), 'o-');
title("Nominal vs. Measured β IDM Median");
xlabel("Nominal Angle (Degrees)");
ylabel("Measured Angle (Degrees)");
xlim([-5 5]);
ylim([-5 5]);
legend("Ideal Case", "Real Case");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

%% Plot cross-track velocities in perfect case and with errors from angle mismeasurements.
%% Using ram velocity of v = 7800 m/s

v = 7800;       % ram velocity magnitude [m/s]

% perfect cross-track velocities @ 7800 m/s ram
perfect_x_dir = -v*tand(sort(alpha_nom));
perfect_y_dir = v*tand(sort(beta_nom));

% measured cross-track velocities @ 7800 m/s ram
measured_x_dir_med = -v*tand(sort(alpha_measured_med));
measured_y_dir_med = v*tand(sort(beta_measured_med));
measured_x_dir_max = -v*tand(sort(alpha_maxes));
measured_y_dir_max = v*tand(sort(beta_maxes));
measured_x_dir_min = -v*tand(sort(alpha_mins));
measured_y_dir_min = v*tand(sort(beta_mins));

% plot v_x
figure(5);
plot(sort(alpha_nom), perfect_x_dir, 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
errorbar(sort(alpha_nom), measured_x_dir_med, measured_x_dir_med - measured_x_dir_min, measured_x_dir_max - measured_x_dir_med, 'o-');
xlabel("Nominal α (Degrees)");
ylabel("v_x (m/s)");
title("Median v_x from Nominal and Measured α at 7800 m/s ram");
xlim([-5 5]);
legend("From Nominal α", "From Measured α Medians");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

% plot v_y
figure(6);
plot(sort(beta_nom), perfect_y_dir, 'LineWidth', 2);
set(gca, 'fontsize', 14);
hold on;
errorbar(sort(beta_nom), measured_y_dir_med, measured_y_dir_med - measured_y_dir_min, measured_y_dir_max - measured_y_dir_med, 'o-');
xlabel("Nominal β (Degrees)");
ylabel("v_y (m/s)");
title("Median v_y from Nominal and Measured β at 7800 m/s ram");
xlim([-5 5]);
legend("From Nominal β", "From Measured β Medians");
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
         'GridColor', 'k', 'MinorGridAlpha', 0.4, 'MinorGridLineStyle', '-', ...
         'MinorGridColor', [0.3,0.3,0.3], 'GridAlpha', 0.4, 'LineWidth', 1.2);

%% Calculate average error between measured medians and nominal values (angles and velocities, angles between -5 and 5 degrees)

disp("All errors calculated for angles in range +/- 5 degrees");

% sort alpha and beta
sorted_alpha_nom = sort(alpha_nom);
sorted_alpha_measured = sort(alpha_measured_med);
sorted_beta_nom = sort(beta_nom);
sorted_beta_measured = sort(beta_measured_med);

% filter velocities so we only include those corresponding to +/- 5 degrees
% of angles.
perfect_x_dir = perfect_x_dir(abs(sorted_alpha_nom) <= 5);
perfect_y_dir = perfect_y_dir(abs(sorted_beta_nom) <= 5);
measured_x_dir_med = measured_x_dir_med(abs(sorted_alpha_nom) <= 5);
measured_y_dir_med = measured_y_dir_med(abs(sorted_beta_nom) <= 5);

% alpha errors
avg_alpha_err = mean(sorted_alpha_nom - sorted_alpha_measured);
avg_alpha_mag_err = mean(abs(sorted_alpha_nom) - abs(sorted_alpha_measured));
max_alpha_mag_err = max(abs(sorted_alpha_nom) - abs(sorted_alpha_measured));
disp(strcat("Average alpha error: ", string(avg_alpha_err), " degrees"));
disp(strcat("Average alpha magnitude error: ", string(avg_alpha_mag_err), " degrees"));
disp(strcat("Maximum alpha magnitude error: ", string(max_alpha_mag_err), " degrees"));
fprintf('\n');

% beta errors
avg_beta_err = mean(sorted_beta_nom - sorted_beta_measured);
avg_beta_mag_err = mean(abs(sorted_beta_nom) - abs(sorted_beta_measured));
max_beta_mag_err = max(abs(sorted_beta_nom) - abs(sorted_beta_measured));
disp(strcat("Average beta error: ", string(avg_beta_err), " degrees"));
disp(strcat("Average beta magnitude error: ", string(avg_beta_mag_err), " degrees"));
disp(strcat("Maximum beta magnitude error: ", string(max_beta_mag_err), " degrees"));
fprintf('\n');

% vx errors
avg_vx_err = mean(perfect_x_dir - measured_x_dir_med);
avg_vx_mag_err = mean(abs(perfect_x_dir) - abs(measured_x_dir_med));
max_vx_mag_err = max(abs(perfect_x_dir) - abs(measured_x_dir_med));
disp(strcat("Average vx error: ", string(avg_vx_err), " m/s"));
disp(strcat("Average vx magnitude error: ", string(avg_vx_mag_err), " m/s"));
disp(strcat("Maximum vx magnitude error: ", string(max_vx_mag_err), " m/s"));
fprintf('\n');

% vy errors
avg_vy_err = mean(perfect_y_dir - measured_y_dir_med);
avg_vy_mag_err = mean(abs(perfect_y_dir) - abs(measured_y_dir_med));
max_vy_mag_err = max(abs(perfect_y_dir) - abs(measured_y_dir_med));
disp(strcat("Average vy error: ", string(avg_vy_err), " m/s"));
disp(strcat("Average vy magnitude error: ", string(avg_vy_mag_err), " m/s"));
disp(strcat("Maximum vy magnitude error: ", string(max_vy_mag_err), " m/s"));
fprintf('\n');

% % average error in alpha
% alpha_err_mag = 0;
% alpha_err = 0;
% max_alpha_err = 0;
% sorted_alpha_nom = sort(alpha_nom);
% sorted_alpha_measured = sort(alpha_measured_med);
% sorted_alpha_nom_pm5degrees = sorted_alpha_nom(abs(sorted_alpha_nom) <= 5);
% sorted_alpha_measured_pm5degrees = sorted_alpha_measured(abs(sorted_alpha_nom) <= 5);
% for i = 1:length(sorted_alpha_nom_pm5degrees)
%     err_mag = abs(sorted_alpha_nom_pm5degrees(i)) - abs(sorted_alpha_measured_pm5degrees(i));
%     err = sorted_alpha_nom_pm5degrees(i) - sorted_alpha_measured_pm5degrees(i);
%     if abs(err_mag) > abs(max_alpha_err)
%         max_alpha_err = err_mag;
%     end
%     alpha_err_mag = alpha_err_mag + err_mag;
%     alpha_err = alpha_err + err;
% end
% alpha_err_mag = alpha_err_mag / length(alpha_nom);
% alpha_err = alpha_err / length(alpha_nom);
% % negative errors indicate overestimation (measured > nominal)
% disp(strcat("Average Alpha Magnitude Error: ", string(alpha_err_mag), " degrees"));
% disp(strcat("Average Alpha Error: ", string(alpha_err), " degrees"));
% disp(strcat("Max Alpha Magnitude Error: ", string(max_alpha_err), " degrees"));
% fprintf('\n');
% 
% % average error in beta
% beta_err_mag = 0;
% beta_err = 0;
% max_beta_err = 0;
% sorted_beta_nom = sort(beta_nom);
% sorted_beta_measured = sort(beta_measured_med);
% sorted_beta_nom_pm5degrees = sorted_beta_nom(abs(sorted_beta_nom) <= 5);
% sorted_beta_measured_pm5degrees = sorted_beta_measured(abs(sorted_beta_nom) <= 5);
% for i = 1:length(sorted_beta_nom_pm5degrees)
%     err_mag = abs(sorted_beta_nom_pm5degrees(i)) - abs(sorted_beta_measured_pm5degrees(i));
%     err = sorted_beta_nom_pm5degrees(i) - sorted_beta_measured_pm5degrees(i);
%     if abs(err_mag) > abs(max_beta_err)
%         max_beta_err = err_mag;
%     end
%     beta_err_mag = beta_err_mag + err_mag;
%     beta_err = beta_err + err;
% end
% beta_err_mag = beta_err_mag / length(beta_nom);
% beta_err = beta_err / length(beta_nom);
% % negative errors indicate overestimation (measured > nominal)
% disp(strcat("Average Beta Magnitude Error: ", string(beta_err_mag), " degrees"));
% disp(strcat("Average Beta Error: ", string(beta_err), " degrees"));
% disp(strcat("Max Beta Magnitude Error: ", string(max_beta_err), " degrees"));
% fprintf('\n');
% 
% % average error in v_x
% vx_err_mag = 0;
% vx_err = 0;
% max_vx_err = 0;
% vx_nom = perfect_x_dir;
% vx_measure = measured_x_dir_med;
% sorted_vx_nom_pm5degrees = vx_nom(abs(sorted_alpha_nom) <= 5);
% sorted_vx_measured_pm5degrees = vx_measure(abs(sorted_alpha_nom) <= 5);
% for i = 1:length(sorted_vx_nom_pm5degrees)
%     err_mag = abs(sorted_vx_nom_pm5degrees(i)) - abs(sorted_vx_measured_pm5degrees(i));
%     err = sorted_vx_nom_pm5degrees(i) - sorted_vx_measured_pm5degrees(i);
%     if abs(err_mag) > abs(max_vx_err)
%         max_vx_err = err_mag;
%     end
%     vx_err_mag = vx_err_mag + err_mag;
%     vx_err = vx_err + err;
% end
% vx_err_mag = vx_err_mag / length(vx_nom);
% vx_err = vx_err / length(vx_nom);
% % negative errors indicate overestimation (measured > nominal)
% disp(strcat("Average v_x Magnitude Error: ", string(vx_err_mag), " m/s"));
% disp(strcat("Average v_x Error: ", string(vx_err), " m/s"));
% disp(strcat("Max v_x Magnitude Error: ", string(max_vx_err), " m/s"));
% disp(alpha_nom())
% fprintf('\n');
% 
% % average error in v_y
% vy_err_mag = 0;
% vy_err = 0;
% max_vy_err = 0;
% vy_nom = perfect_y_dir;
% vy_measure = measured_y_dir_med;
% sorted_vy_nom_pm5degrees = vy_nom(abs(sorted_beta_nom) <= 5);
% sorted_vy_measured_pm5degrees = vy_measure(abs(sorted_beta_nom) <= 5);
% for i = 1:length(sorted_vy_nom_pm5degrees)
%     err_mag = abs(sorted_vy_nom_pm5degrees(i)) - abs(sorted_vy_measured_pm5degrees(i));
%     err = sorted_vy_nom_pm5degrees(i) - sorted_vy_measured_pm5degrees(i);
%     if abs(err_mag) > abs(max_vy_err)
%         max_vy_err = err_mag;
%     end
%     vy_err_mag = vy_err_mag + err_mag;
%     vy_err = vy_err + err;
% end
% vy_err_mag = vy_err_mag / length(vy_nom);
% vy_err = vy_err / length(vy_nom);
% % negative errors indicate overestimation (measured > nominal)
% disp(strcat("Average v_y Magnitude Error: ", string(vy_err_mag), " m/s"));
% disp(strcat("Average v_y Error: ", string(vy_err), " m/s"));
% disp(strcat("Max v_y Magnitude Error: ", string(max_vy_err), " m/s"));
% fprintf('\n');
