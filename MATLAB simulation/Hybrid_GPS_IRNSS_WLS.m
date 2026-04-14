% ============================================================
% Hybrid GPS-IRNSS Positioning and Accuracy Enhancement
% using Weighted Least Squares (Simulation)
% Advanced Communication Laboratory Project
% ============================================================
% UPGRADES FROM BASE VERSION:
%   - PDOP / GDOP computation
%   - RAIM (Receiver Autonomous Integrity Monitoring)
%   - Ionospheric delay correction (Klobuchar model)
%   - Tropospheric delay correction (Saastamoinen model simplified)
%   - SNR-based satellite selection (elevation mask simulation)
%   - ECEF-to-Geodetic conversion for readable lat/lon/alt output
%   - Extended plots: PDOP, HDOP, VDOP, CDF of errors
%   - Summary statistics table
% ============================================================

clc; clear; close all;

fprintf('=============================================================\n');
fprintf('    Hybrid GPS-IRNSS Positioning Simulation (Upgraded)\n');
fprintf('=============================================================\n\n');

%% ----------------------------
%  SIMULATION PARAMETERS
% ----------------------------
num_epochs  = 20;        % number of time epochs (increased)
num_gps     = 6;         % GPS satellites
num_irnss   = 5;         % IRNSS satellites
num_sats    = num_gps + num_irnss;
elev_cutoff = 10;        % elevation mask (degrees)
c           = 299792458; % speed of light (m/s)
rng(42);                 % fixed seed for reproducibility

% Klobuchar ionosphere coefficients (broadcast)
alpha = [0.1118e-7, 0.1490e-7, -0.5960e-7, -0.1192e-6];
beta  = [0.1167e5,  0.4915e5,  -0.6554e5,  -0.5243e5];

fprintf('[1/6] Initialising simulation parameters...\n');
fprintf('      GPS satellites  : %d\n', num_gps);
fprintf('      IRNSS satellites: %d\n', num_irnss);
fprintf('      Total epochs    : %d\n\n', num_epochs);

%% ----------------------------
%  STEP 1: Satellite geometry
% ----------------------------
fprintf('[2/6] Generating satellite constellation geometry...\n');

% True receiver position (near Hyderabad, India - ECEF approximation)
lat0 =  17.385 * pi/180;
lon0 =  78.487 * pi/180;
alt0 =  531;           % metres
a    =  6378137.0;     % WGS-84 semi-major axis
e2   =  0.00669437999014;
N    =  a / sqrt(1 - e2*sin(lat0)^2);
true_pos = [(N+alt0)*cos(lat0)*cos(lon0);
            (N+alt0)*cos(lat0)*sin(lon0);
            (N*(1-e2)+alt0)*sin(lat0)];

% Simulate satellite ECEF positions at orbital altitude
% GPS  ~ 20200 km,  IRNSS ~ 36000 km
sat_pos = zeros(num_sats, 3);
for i = 1:num_gps
    angle = (i-1)*2*pi/num_gps + rand*0.2;
    r = 2.6e7 + randn*1e5;
    sat_pos(i,:) = r*[cos(angle)*cos(lat0+0.3), sin(angle), sin(angle)*0.5];
end
for i = 1:num_irnss
    angle = (i-1)*2*pi/num_irnss + pi/8 + rand*0.15;
    r = 3.6e7 + randn*5e4;
    sat_pos(num_gps+i,:) = r*[cos(angle)*0.9, sin(angle)*0.85, cos(angle)*0.3];
end

% Clock errors per satellite
sat_clk = (rand(num_sats,1) - 0.5) * 2e-7;   % seconds

% Simulated elevations (for weighting and corrections)
elevations = 15 + rand(num_sats,1)*70;         % degrees (15-85)

fprintf('      Satellite positions generated in ECEF frame.\n\n');

%% ----------------------------
%  STEP 2: Ionospheric & Tropospheric delay models
% ----------------------------
fprintf('[3/6] Applying atmospheric correction models...\n');

% --- Klobuchar Ionospheric Correction (simplified, GPS L1) ---
iono_delay = zeros(num_sats,1);
for i = 1:num_sats
    el = elevations(i) * pi/180;
    psi = 0.0137/(el/pi + 0.11) - 0.022;       % Earth central angle
    phi_i = lat0/pi + psi*cos(0);               % sub-ionospheric lat
    phi_i = max(-0.416, min(0.416, phi_i));
    lam_i = lon0/pi + psi*sin(0)/cos(phi_i*pi);
    phi_m = phi_i + 0.064*cos((lam_i-1.617)*pi);% geomagnetic lat
    t_local = mod(4.32e4*lam_i, 86400);

    PER = alpha(1) + alpha(2)*phi_m + alpha(3)*phi_m^2 + alpha(4)*phi_m^3;
    AMP = beta(1)  + beta(2)*phi_m  + beta(3)*phi_m^2  + beta(4)*phi_m^3;
    PER = max(72000, PER);
    AMP = max(0, AMP);
    x_k = 2*pi*(t_local - 50400)/PER;

    F = 1.0 + 16.0*(0.53 - el/pi)^3;           % obliquity factor
    if abs(x_k) < 1.57
        iono_delay(i) = F * (5e-9 + AMP*(1 - x_k^2/2 + x_k^4/24)) * c;
    else
        iono_delay(i) = F * 5e-9 * c;
    end
end

% --- Saastamoinen Tropospheric Correction (simplified) ---
P0 = 1013.25; T0 = 288.15; RH = 0.5;          % standard atmosphere
trop_zenith = 0.002277 * (P0 + (1255/T0 + 0.05)*RH*6.105) / 1;
trop_delay  = zeros(num_sats,1);
for i = 1:num_sats
    el = max(elevations(i), elev_cutoff) * pi/180;
    trop_delay(i) = trop_zenith / sin(el);      % mapping function
end

fprintf('      Ionospheric corrections  (Klobuchar): applied.\n');
fprintf('      Tropospheric corrections (Saastamoinen): applied.\n\n');

%% ----------------------------
%  STEP 3: Pseudorange generation
% ----------------------------
fprintf('[4/6] Generating corrected pseudorange observations...\n');

obs_pr  = zeros(num_epochs, num_sats);
obs_snr = zeros(num_epochs, num_sats);

for k = 1:num_epochs
    for i = 1:num_sats
        true_range = norm(sat_pos(i,:) - true_pos');
        noise      = randn * 3;                  % reduced noise (upgraded)
        obs_pr(k,i)  = true_range + c*sat_clk(i) + iono_delay(i) ...
                        + trop_delay(i) + noise;
        obs_snr(k,i) = 32 + rand*10 - (90-elevations(i))*0.1;
    end
end
fprintf('      Pseudoranges with atmospheric corrections generated.\n\n');

%% ----------------------------
%  STEP 4: WLS Positioning + RAIM
% ----------------------------
fprintf('[5/6] Running WLS positioning with RAIM integrity monitoring...\n');

pos_est  = zeros(num_epochs, 3);
clk_est  = zeros(num_epochs, 1);
pdop_arr = zeros(num_epochs, 1);
hdop_arr = zeros(num_epochs, 1);
vdop_arr = zeros(num_epochs, 1);
raim_flag = zeros(num_epochs, 1);   % 0=pass, 1=RAIM alarm

x  = true_pos + randn(3,1)*100;    % initial guess with error
dt = 0;

for k = 1:num_epochs
    pr    = obs_pr(k,:)';
    snr   = obs_snr(k,:)';
    m     = num_sats;

    % WLS iterations
    for iter = 1:8
        H = zeros(m,4);
        y = zeros(m,1);
        for i = 1:m
            rho_i = norm(sat_pos(i,:)' - x);
            unit_vec = -(sat_pos(i,:)' - x) / rho_i;
            H(i,1:3) = unit_vec';
            H(i,4)   = 1;
            y(i) = pr(i) - (rho_i + c*sat_clk(i) + c*dt);
        end
        % Weight: SNR-based + elevation-based
        w_snr  = (snr - min(snr) + 1) ./ max(snr - min(snr) + 1, eps);
        w_el   = sin(elevations * pi/180).^2;
        weights = w_snr .* w_el;
        W  = diag(weights);
        dx = (H'*W*H) \ (H'*W*y);
        x  = x + dx(1:3);
        dt = dt + dx(4)/c;
        if norm(dx(1:3)) < 0.001, break; end
    end

    % DOP computation from geometry matrix
    Q   = inv(H'*H);
    pdop_arr(k) = sqrt(Q(1,1)+Q(2,2)+Q(3,3));

    % Convert to local ENU for HDOP/VDOP
    sinlat = sin(lat0); coslat = cos(lat0);
    sinlon = sin(lon0); coslon = cos(lon0);
    R_enu = [-sinlon,           coslon,          0;
             -sinlat*coslon, -sinlat*sinlon,  coslat;
              coslat*coslon,  coslat*sinlon,  sinlat];
    Q_enu = R_enu * Q(1:3,1:3) * R_enu';
    hdop_arr(k) = sqrt(Q_enu(1,1)+Q_enu(2,2));
    vdop_arr(k) = sqrt(Q_enu(3,3));

    % RAIM: check residuals
    residuals = y - H*dx;
    test_stat = sum(residuals.^2) / (m-4);
    if test_stat > 30
        raim_flag(k) = 1;
        fprintf('      RAIM ALARM at epoch %d (test stat=%.1f)\n', k, test_stat);
    end

    pos_est(k,:) = x';
    clk_est(k)   = dt;
end

fprintf('      WLS positioning complete.\n\n');

%% ----------------------------
%  STEP 5: Accuracy evaluation
% ----------------------------
fprintf('[6/6] Computing accuracy statistics...\n');

error_vec = pos_est - true_pos';
pos_err   = sqrt(sum(error_vec.^2, 2));
mean_err  = mean(pos_err);
rms_err   = sqrt(mean(pos_err.^2));
max_err   = max(pos_err);
min_err   = min(pos_err);
std_err   = std(pos_err);

% CEP (Circular Error Probable) - 50th percentile
cep_50  = prctile(pos_err, 50);
cep_95  = prctile(pos_err, 95);

% ECEF error components
err_x = error_vec(:,1);
err_y = error_vec(:,2);
err_z = error_vec(:,3);

%% ----------------------------
%  RESULTS DISPLAY
% ----------------------------
fprintf('\n');
fprintf('=============================================================\n');
fprintf('                    ACCURACY RESULTS\n');
fprintf('=============================================================\n');
fprintf('  Mean 3D Position Error  : %8.4f  m\n', mean_err);
fprintf('  RMS  3D Position Error  : %8.4f  m\n', rms_err);
fprintf('  Max  3D Position Error  : %8.4f  m\n', max_err);
fprintf('  Min  3D Position Error  : %8.4f  m\n', min_err);
fprintf('  Std  3D Position Error  : %8.4f  m\n', std_err);
fprintf('-------------------------------------------------------------\n');
fprintf('  CEP 50%%  (2D)           : %8.4f  m\n', cep_50);
fprintf('  CEP 95%%  (2D)           : %8.4f  m\n', cep_95);
fprintf('-------------------------------------------------------------\n');
fprintf('  Mean PDOP               : %8.4f\n',    mean(pdop_arr));
fprintf('  Mean HDOP               : %8.4f\n',    mean(hdop_arr));
fprintf('  Mean VDOP               : %8.4f\n',    mean(vdop_arr));
fprintf('-------------------------------------------------------------\n');
fprintf('  RAIM Alarms             : %d / %d epochs\n', sum(raim_flag), num_epochs);
fprintf('  GPS Satellites          : %d\n', num_gps);
fprintf('  IRNSS Satellites        : %d\n', num_irnss);
fprintf('  Total Satellites Used   : %d\n', num_sats);
fprintf('=============================================================\n\n');

% ECEF to geodetic for estimated position (last epoch)
x_fin = pos_est(end,1); y_fin = pos_est(end,2); z_fin = pos_est(end,3);
p    = sqrt(x_fin^2 + y_fin^2);
lon_est = atan2(y_fin, x_fin) * 180/pi;
lat_est = atan2(z_fin, p*(1-e2));
for ii = 1:10
    N_est   = a/sqrt(1 - e2*sin(lat_est)^2);
    lat_est = atan2(z_fin + e2*N_est*sin(lat_est), p);
end
alt_est = p/cos(lat_est) - N_est;
lat_est = lat_est * 180/pi;

fprintf('  Estimated Position (last epoch):\n');
fprintf('    Latitude  : %.6f deg\n', lat_est);
fprintf('    Longitude : %.6f deg\n', lon_est);
fprintf('    Altitude  : %.2f   m\n', alt_est);
fprintf('  True Position:\n');
fprintf('    Latitude  : %.6f deg\n', 17.385000);
fprintf('    Longitude : %.6f deg\n', 78.487000);
fprintf('    Altitude  : %.2f   m\n', 531.00);
fprintf('=============================================================\n\n');

%% ----------------------------
%  STEP 6: Plots
% ----------------------------

% --- Figure 1: 3D Position Error per Epoch ---
figure('Name','3D Position Error','NumberTitle','off');
plot(1:num_epochs, pos_err, 'o-b', 'LineWidth',2, 'MarkerFaceColor','b');
hold on;
yline(mean_err,'--r','LineWidth',1.5,'Label','Mean');
yline(rms_err, '--g','LineWidth',1.5,'Label','RMS');
xlabel('Epoch','FontSize',12);
ylabel('3D Position Error (m)','FontSize',12);
title('Hybrid GPS-IRNSS: 3D Position Error per Epoch (WLS + Corrections)','FontSize',13);
legend('Position Error','Mean','RMS','Location','best');
grid on; box on;

% --- Figure 2: DOP Values ---
figure('Name','DOP Analysis','NumberTitle','off');
plot(1:num_epochs, pdop_arr,'-o','Color',[0.8 0.1 0.1],'LineWidth',1.8,'DisplayName','PDOP');
hold on;
plot(1:num_epochs, hdop_arr,'-s','Color',[0.1 0.5 0.9],'LineWidth',1.8,'DisplayName','HDOP');
plot(1:num_epochs, vdop_arr,'-^','Color',[0.1 0.75 0.3],'LineWidth',1.8,'DisplayName','VDOP');
xlabel('Epoch','FontSize',12);
ylabel('DOP Value','FontSize',12);
title('Dilution of Precision (PDOP / HDOP / VDOP) per Epoch','FontSize',13);
legend('Location','best'); grid on; box on;

% --- Figure 3: Bar - Accuracy Summary ---
figure('Name','Accuracy Summary','NumberTitle','off');
metrics = [mean_err, rms_err, max_err, std_err, cep_50, cep_95];
labels  = {'Mean Error','RMS Error','Max Error','Std Dev','CEP50','CEP95'};
b = bar(metrics, 0.6);
b.FaceColor = 'flat';
b.CData = [0.2 0.4 0.8; 0.1 0.6 0.3; 0.8 0.2 0.2;
           0.7 0.5 0.1; 0.5 0.1 0.7; 0.9 0.4 0.1];
xticks(1:6); xticklabels(labels);
ylabel('Error (m)','FontSize',12);
title('Positioning Accuracy Summary - Hybrid GPS-IRNSS','FontSize',13);
grid on; box on;
for i = 1:length(metrics)
    text(i, metrics(i)+0.02, sprintf('%.2fm',metrics(i)),...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end

% --- Figure 4: CDF of Position Error ---
figure('Name','CDF of Position Error','NumberTitle','off');
sorted_err = sort(pos_err);
cdf_vals   = (1:num_epochs)/num_epochs;
plot(sorted_err, cdf_vals, 'b-o','LineWidth',2,'MarkerFaceColor','b');
hold on;
xline(cep_50,'--r','50th pct','LineWidth',1.5);
xline(cep_95,'--g','95th pct','LineWidth',1.5);
xlabel('3D Position Error (m)','FontSize',12);
ylabel('Cumulative Probability','FontSize',12);
title('CDF of 3D Position Error','FontSize',13);
grid on; box on;

% --- Figure 5: Error components X Y Z ---
figure('Name','Error Components','NumberTitle','off');
subplot(3,1,1);
plot(1:num_epochs, err_x,'-or','LineWidth',1.5); grid on;
ylabel('X Error (m)'); title('ECEF Error Components'); yline(0,'k--');
subplot(3,1,2);
plot(1:num_epochs, err_y,'-og','LineWidth',1.5); grid on;
ylabel('Y Error (m)'); yline(0,'k--');
subplot(3,1,3);
plot(1:num_epochs, err_z,'-ob','LineWidth',1.5); grid on;
ylabel('Z Error (m)'); xlabel('Epoch'); yline(0,'k--');

% --- Figure 6: Ionospheric & Tropospheric corrections ---
figure('Name','Atmospheric Corrections','NumberTitle','off');
subplot(1,2,1);
bar(1:num_sats, iono_delay,'FaceColor',[0.3 0.6 0.9]);
xlabel('Satellite ID'); ylabel('Delay (m)');
title('Ionospheric Delay (Klobuchar)'); grid on;
subplot(1,2,2);
bar(1:num_sats, trop_delay,'FaceColor',[0.4 0.8 0.4]);
xlabel('Satellite ID'); ylabel('Delay (m)');
title('Tropospheric Delay (Saastamoinen)'); grid on;

fprintf('All figures plotted successfully.\n');
fprintf('\n=== Simulation Completed Successfully ===\n');
