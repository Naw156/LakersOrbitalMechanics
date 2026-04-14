clc; clear; close all

%% Constants
muE = 398600.4418;      % km^3/s^2
RE  = 6378.1363;        % km
global mu 
mu = muE;

%% ---------------- PART 1: DEFINE PARKING ORBIT ----------------
h_park   = 185;         % km altitude
inc_deg  = 28.5;        % deg
RAAN_deg = 0;           % deg
w_deg    = 0;           % deg
theta    = 0;           % deg

rpark = RE + h_park;    % km
epark = 0;              % circular orbit
hmag  = sqrt(muE * rpark);

% Initial Cartesian state in geocentric equatorial frame
[r0_vec, v0_vec] = Perifocal2GE(hmag, inc_deg, RAAN_deg, epark, w_deg, theta, muE);

fprintf('Initial parking-orbit state:\n')
fprintf('r0 = [%.6f %.6f %.6f] km\n', r0_vec(1), r0_vec(2), r0_vec(3))
fprintf('v0 = [%.6f %.6f %.6f] km/s\n\n', v0_vec(1), v0_vec(2), v0_vec(3))

% Recover orbital elements as a check
[h0,i0,W0,e0,w0,th0] = OrbitalElements(r0_vec,v0_vec,muE);

% Semi-major axis from h and e
p0 = h0^2 / muE;
a0 = p0 / (1 - e0^2);

% Orbital period
Torb = 2*pi*sqrt(a0^3/muE);

fprintf('Recovered orbital elements:\n')
fprintf('h     = %.6f km^2/s\n', h0)
fprintf('i     = %.6f deg\n', i0)
fprintf('RAAN  = %.6f deg\n', W0)
fprintf('e     = %.10f\n', e0)
fprintf('w     = %.6f deg\n', w0)
fprintf('theta = %.6f deg\n\n', th0)

fprintf('Semi-major axis a = %.6f km\n', a0)
fprintf('Orbital period T  = %.6f s (%.6f min)\n\n', Torb, Torb/60)

%% Time grid for one full orbit
N = 500;
tspan = linspace(0, Torb, N);

%% ---------------- PART 2A: ANALYTICAL PROPAGATION ----------------
r_analytical = zeros(N,3);
v_analytical = zeros(N,3);

for k = 1:N
    [rtemp, vtemp] = UVars(r0_vec, v0_vec, tspan(k));
    r_analytical(k,:) = rtemp(:).';
    v_analytical(k,:) = vtemp(:).';
end

%% ---------------- PART 2B: NUMERICAL PROPAGATION ----------------
dydt = @(t,y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3];
y0 = [r0_vec; v0_vec];

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

[t_num, y_num] = ode45(dydt, tspan, y0, opts);

r_numerical = y_num(:,1:3);
v_numerical = y_num(:,4:6);

%% Final-state comparison after one orbit
rA_end = r_analytical(end,:).';
vA_end = v_analytical(end,:).';

rN_end = r_numerical(end,:).';
vN_end = v_numerical(end,:).';

fprintf('--- Final state after one orbit ---\n')
fprintf('Analytical:\n')
fprintf('rA(T) = [%.6f %.6f %.6f] km\n', rA_end(1), rA_end(2), rA_end(3))
fprintf('vA(T) = [%.6f %.6f %.6f] km/s\n\n', vA_end(1), vA_end(2), vA_end(3))

fprintf('Numerical:\n')
fprintf('rN(T) = [%.6f %.6f %.6f] km\n', rN_end(1), rN_end(2), rN_end(3))
fprintf('vN(T) = [%.6f %.6f %.6f] km/s\n\n', vN_end(1), vN_end(2), vN_end(3))

% Error relative to initial state
err_rA = norm(rA_end - r0_vec(:));
err_vA = norm(vA_end - v0_vec(:));

err_rN = norm(rN_end - r0_vec(:));
err_vN = norm(vN_end - v0_vec(:));

% Difference between methods
err_rAN = norm(rA_end - rN_end);
err_vAN = norm(vA_end - vN_end);

fprintf('Error relative to initial state:\n')
fprintf('Analytical: |dr| = %.12e km,    |dv| = %.12e km/s\n', err_rA, err_vA)
fprintf('Numerical : |dr| = %.12e km,    |dv| = %.12e km/s\n\n', err_rN, err_vN)

fprintf('Difference between analytical and numerical final states:\n')
fprintf('|dr| = %.12e km\n', err_rAN)
fprintf('|dv| = %.12e km/s\n', err_vAN)

%% Plots
figure
plot3(r_analytical(:,1), r_analytical(:,2), r_analytical(:,3), 'b-', 'LineWidth', 1.5)
hold on
plot3(r_numerical(:,1), r_numerical(:,2), r_numerical(:,3), 'r--', 'LineWidth', 1.2)
plot3(r0_vec(1), r0_vec(2), r0_vec(3), 'ko', 'MarkerFaceColor', 'k')

grid on
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('Analytical (UVars)', 'Numerical (ode45)', 'Initial state', 'Earth', 'Location', 'best')
title('Parking Orbit Propagation for One Full Revolution')
view(3)

figure
plot(tspan/60, vecnorm(r_analytical,2,2), 'b-', 'LineWidth', 1.5)
hold on
plot(t_num/60, vecnorm(r_numerical,2,2), 'r--', 'LineWidth', 1.2)
grid on
xlabel('Time (min)')
ylabel('|r| (km)')
legend('Analytical', 'Numerical', 'Location', 'best')
title('Radius Magnitude vs Time')

%%%%%%%%%% TLI
%% ---------------- PART 3: HOHMANN-LIKE TRANSLUNAR INJECTION ----------------
muM = 4902.800066;      % km^3/s^2, Moon
RM  = 1737.4;           % km

% Mean Earth-Moon distance for first-order patched-conic estimate
r_moon = 384400;        % km

% Initial orbit: circular parking orbit
rA  = rpark;
rAp = rpark;

% Target orbit for Hohmann function:
% circular Earth-centered orbit at Moon distance
rB  = r_moon;
rBp = r_moon;

[delV1, delV2, t_trans] = Hohmann(rA, rAp, rB, rBp, muE);

fprintf('\n--- Hohmann-like Translunar Injection ---\n')
fprintf('Parking orbit radius           = %.6f km\n', rpark)
fprintf('Target radius                  = %.6f km\n', r_moon)
fprintf('TLI Delta-V                    = %.6f km/s\n', delV1)
fprintf('Hohmann circularization Delta-V= %.6f km/s\n', delV2)
fprintf('Transfer time to apogee        = %.6f s\n', t_trans)
fprintf('Transfer time to apogee        = %.6f hr\n', t_trans/3600)
fprintf('Transfer time to apogee        = %.6f days\n', t_trans/86400)


%% ---------------- PART 8: LAUNCH-TIME PHASING SWEEP ----------------
% Sweep launch time forward only from the nominal date
dt_hours = 0 : 0.1 : 8;
nSweep   = length(dt_hours);

miss_distance = zeros(nSweep,1);
v_inf_sweep   = zeros(nSweep,1);
    launch_date0 = datetime(2026,4,1,12,0,0);
launch_dates  = launch_date0 + hours(dt_hours);

r_moon_arr_all = zeros(nSweep,3);
r_sc_arr_all   = zeros(nSweep,3);

for j = 1:nSweep
    dt_sec = dt_hours(j) * 3600;
    
    [~, ~, ~, ~, r_sc_arr_j, v_sc_arr_j, r_moon_arr_j, v_moon_arr_j, ...
     miss_distance(j), v_inf_sweep(j)] = ...
        PropagateTLI(r0_vec, v0_vec, launch_date0, dt_sec, delV1, t_trans, ...
                     muE, muM, 2);
    
    r_sc_arr_all(j,:)   = r_sc_arr_j.';
    r_moon_arr_all(j,:) = r_moon_arr_j.';
end

[miss_best, idx_best] = min(miss_distance);
best_launch_date  = launch_dates(idx_best);
best_arrival_date = best_launch_date + seconds(t_trans);

fprintf('\n--- Launch-Time Phasing Sweep ---\n')
fprintf('Nominal launch date           = %s\n', datestr(launch_date0))
fprintf('Best launch date in sweep     = %s\n', datestr(best_launch_date))
fprintf('Best arrival date             = %s\n', datestr(best_arrival_date))
fprintf('Minimum miss distance         = %.6f km\n', miss_best)
fprintf('Arrival v_infinity there      = %.6f km/s\n', v_inf_sweep(idx_best))

%% ---------------- PART 9: REPORT BEST-CASE ARRIVAL GEOMETRY ----------------
r_sc_best   = r_sc_arr_all(idx_best,:).';
r_moon_best = r_moon_arr_all(idx_best,:).';

fprintf('\nBest-case arrival vectors:\n')
fprintf('Spacecraft arrival position   = [%.6f %.6f %.6f] km\n', ...
    r_sc_best(1), r_sc_best(2), r_sc_best(3))
fprintf('Moon arrival position         = [%.6f %.6f %.6f] km\n', ...
    r_moon_best(1), r_moon_best(2), r_moon_best(3))
fprintf('Best miss distance            = %.6f km\n', norm(r_sc_best - r_moon_best))
%% ---------------- PART 9B: PROPAGATE AND PLOT OPTIMIZED TRAJECTORY ----------------
N_opt = 1200;
dt_best_sec = dt_hours(idx_best) * 3600;

[r_depart_best, v_depart_best, r_hist_best, v_hist_best, ...
 r_arr_best, v_arr_best, r_moon_arr_best, v_moon_arr_best, ...
 miss_arrival_best, v_inf_best, t_hist_best] = ...
    PropagateTLI(r0_vec, v0_vec, launch_date0, dt_best_sec, delV1, t_trans, ...
                 muE, muM, N_opt);

% Basic optimized-burn properties
v_park_best = norm(v_depart_best - delV1 * (v_depart_best / norm(v_depart_best)));
v_depart_mag = norm(v_depart_best);
C3_best = v_depart_mag^2 - 2*muE/norm(r_depart_best);

fprintf('\n--- Optimized TLI Properties ---\n')
fprintf('Best launch date              = %s\n', datestr(best_launch_date))
fprintf('Best arrival date             = %s\n', datestr(best_arrival_date))
fprintf('Launch delay from nominal     = %.6f hr\n', dt_hours(idx_best))
fprintf('TLI Delta-V                   = %.6f km/s\n', delV1)
fprintf('Departure radius              = %.6f km\n', norm(r_depart_best))
fprintf('Departure speed after TLI     = %.6f km/s\n', v_depart_mag)
fprintf('Arrival v_infinity            = %.6f km/s\n', v_inf_best)
fprintf('Arrival miss distance(center) = %.6f km\n', miss_arrival_best)
fprintf('Arrival miss distance(surface)= %.6f km\n', miss_arrival_best - RM)
fprintf('Approx characteristic energy C3 = %.6f km^2/s^2\n', C3_best)

%% Compute true closest approach to the Moon over the optimized trajectory
nHist = length(t_hist_best);
r_moon_hist_best = zeros(nHist,3);
moon_sep_best = zeros(nHist,1);

jd_launch_best = juliandate(best_launch_date);

for k = 1:nHist
    jd_k = jd_launch_best + t_hist_best(k)/86400;
    [rMoon_k, ~] = planetEphemeris(jd_k,'Earth','Moon');
    rMoon_k = rMoon_k(:);

    r_moon_hist_best(k,:) = rMoon_k.';
    moon_sep_best(k) = norm(r_hist_best(k,:).' - rMoon_k);
end

[min_sep_center, idx_ca] = min(moon_sep_best);
min_sep_surface = min_sep_center - RM;
t_ca = t_hist_best(idx_ca);
r_sc_ca = r_hist_best(idx_ca,:).';
r_moon_ca = r_moon_hist_best(idx_ca,:).';

fprintf('\n--- Closest Approach Along Optimized Trajectory ---\n')
fprintf('Closest approach time after TLI = %.6f hr\n', t_ca/3600)
fprintf('Closest distance to Moon center = %.6f km\n', min_sep_center)
fprintf('Closest altitude above surface  = %.6f km\n', min_sep_surface)
fprintf('Spacecraft position at CA       = [%.6f %.6f %.6f] km\n', ...
    r_sc_ca(1), r_sc_ca(2), r_sc_ca(3))
fprintf('Moon position at CA             = [%.6f %.6f %.6f] km\n', ...
    r_moon_ca(1), r_moon_ca(2), r_moon_ca(3))

%% Plot optimized trajectory
figure
plot3(r_hist_best(:,1), r_hist_best(:,2), r_hist_best(:,3), 'b-', 'LineWidth', 1.5)
hold on
plot3(r_moon_hist_best(:,1), r_moon_hist_best(:,2), r_moon_hist_best(:,3), 'm--', 'LineWidth', 1.2)
plot3(r_depart_best(1), r_depart_best(2), r_depart_best(3), 'go', 'MarkerFaceColor', 'g')
plot3(r_arr_best(1), r_arr_best(2), r_arr_best(3), 'ko', 'MarkerFaceColor', 'k')
plot3(r_sc_ca(1), r_sc_ca(2), r_sc_ca(3), 'ro', 'MarkerFaceColor', 'r')
plot3(r_moon_ca(1), r_moon_ca(2), r_moon_ca(3), 'mo', 'MarkerFaceColor', 'm')

[xE,yE,zE] = sphere(100);
surf(RE*xE, RE*yE, RE*zE, ...
    'FaceColor', [0.6 0.8 1.0], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.25)

grid on
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('Optimized spacecraft trajectory', 'Moon trajectory', ...
       'TLI point', 'Arrival point', 'Closest approach spacecraft point', ...
       'Closest approach Moon point', 'Earth', 'Location', 'best')
title('Optimized Translunar Trajectory')
view(3)

figure
plot(t_hist_best/3600, moon_sep_best, 'b-', 'LineWidth', 1.5)
hold on
yline(RM, 'k--', 'Moon radius')
grid on
xlabel('Time after TLI (hr)')
ylabel('Distance to Moon center (km)')
title('Moon-Relative Separation Along Optimized Trajectory')

figure
plot(t_hist_best/3600, moon_sep_best - RM, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Time after TLI (hr)')
ylabel('Altitude above lunar surface (km)')
title('Altitude Above Lunar Surface Along Optimized Trajectory')

%% ---------------- PART 10: PLOTS FOR PHASING SWEEP ----------------
figure
plot(dt_hours, miss_distance, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Launch time offset from nominal (hr)')
ylabel('Miss distance at arrival (km)')
title('Moon Miss Distance vs Launch-Time Offset')

figure
plot(dt_hours, v_inf_sweep, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Launch time offset from nominal (hr)')
ylabel('Arrival v_{\infty} (km/s)')
title('Arrival v_{\infty} vs Launch-Time Offset')

figure
plot3(r_sc_arr_all(:,1), r_sc_arr_all(:,2), r_sc_arr_all(:,3), 'b.', 'MarkerSize', 10)
hold on
plot3(r_moon_arr_all(:,1), r_moon_arr_all(:,2), r_moon_arr_all(:,3), 'r.', 'MarkerSize', 10)
plot3(r_sc_best(1), r_sc_best(2), r_sc_best(3), 'ko', 'MarkerFaceColor', 'k')
plot3(r_moon_best(1), r_moon_best(2), r_moon_best(3), 'mo', 'MarkerFaceColor', 'm')

[xE,yE,zE] = sphere(100);
surf(RE*xE, RE*yE, RE*zE, ...
    'FaceColor', [0.6 0.8 1.0], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.20)

grid on
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('Spacecraft arrival points','Moon positions at arrival', ...
       'Best spacecraft point','Best Moon point','Earth', ...
       'Location','best')
title('Arrival Geometry During Launch-Time Sweep')
view(3)

%% ---------------- LOCAL FUNCTIONS ----------------
function [r_depart, v_depart, r_hist, v_hist, r_arr, v_arr, ...
          r_moon_arr, v_moon_arr, miss_distance, v_inf_arr, t_hist] = ...
    PropagateTLI(r0_vec, v0_vec, launch_date0, dt_sec, delV1, t_trans, ...
                 muE, muM, Nout)

    % Propagate parking orbit to delayed TLI time using two-body motion
    [r_park, v_park] = UVars(r0_vec, v0_vec, dt_sec);
    r_park = r_park(:);
    v_park = v_park(:);

    % Apply tangential TLI burn
    vhat = v_park / norm(v_park);
    r_depart = r_park;
    v_depart = v_park + delV1 * vhat;

    % Earth + Moon gravity propagation from TLI to arrival
    launch_date = launch_date0 + seconds(dt_sec);
    jd0 = juliandate(launch_date);

    y0 = [r_depart; v_depart];
    tspan = linspace(0, t_trans, Nout);
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

    [t_hist, y_hist] = ode45(@(t,y) stateMoonGravity(t,y,muE,muM,jd0), ...
                             tspan, y0, opts);

    r_hist = y_hist(:,1:3);
    v_hist = y_hist(:,4:6);

    r_arr = r_hist(end,:).';
    v_arr = v_hist(end,:).';

    % Moon state at arrival
    arrival_date = launch_date + seconds(t_trans);
    jd_arrival = juliandate(arrival_date);
    [r_moon_arr, v_moon_arr] = planetEphemeris(jd_arrival,'Earth','Moon');
    r_moon_arr = r_moon_arr(:);
    v_moon_arr = v_moon_arr(:);

    % Metrics
    miss_distance = norm(r_arr - r_moon_arr);
    v_inf_arr     = norm(v_arr - v_moon_arr);
end

function dydt = stateMoonGravity(t,y,muE,muM,jd0)

    r = y(1:3);
    v = y(4:6);

    jd = jd0 + t/86400;

    [rMoon, ~] = planetEphemeris(jd,'Earth','Moon');
    rMoon = rMoon(:);

    Re = r;
    Rm = r - rMoon;

    a = -muE * Re / norm(Re)^3 ...
        -muM * Rm / norm(Rm)^3;

    dydt = [v; a];
end
