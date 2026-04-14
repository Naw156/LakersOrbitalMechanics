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
theta  = 0;           % deg

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

% Earth
%[xE,yE,zE] = sphere(100);
%surf(RE*xE, RE*yE, RE*zE, ...
 %   'FaceColor', [0.6 0.8 1.0], ...
   % 'EdgeColor', 'none', ...
   % 'FaceAlpha', 0.8)

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

%% ---------------- PART 4: POST-TLI STATE ----------------
% Apply the TLI burn tangentially to the current parking-orbit velocity
vhat0 = v0_vec / norm(v0_vec);

rTLI_vec = r0_vec;
vTLI_vec = v0_vec + delV1 * vhat0;

fprintf('\nPost-TLI departure state:\n')
fprintf('rTLI = [%.6f %.6f %.6f] km\n', rTLI_vec(1), rTLI_vec(2), rTLI_vec(3))
fprintf('vTLI = [%.6f %.6f %.6f] km/s\n', vTLI_vec(1), vTLI_vec(2), vTLI_vec(3))

%% ---------------- PART 5: TRANSFER ORBIT CHECKS ----------------
a_trans = (rpark + r_moon)/2;
e_trans = (r_moon - rpark)/(r_moon + rpark);

v_circ_LEO      = sqrt(muE/rpark);
v_perigee_trans = sqrt(muE*(2/rpark - 1/a_trans));
v_apogee_trans  = sqrt(muE*(2/r_moon - 1/a_trans));
v_moon_circ     = sqrt(muE/r_moon);

v_inf_moon = abs(v_moon_circ - v_apogee_trans);

fprintf('\nTransfer-orbit properties:\n')
fprintf('a_transfer                    = %.6f km\n', a_trans)
fprintf('e_transfer                    = %.10f\n', e_trans)
fprintf('LEO circular speed            = %.6f km/s\n', v_circ_LEO)
fprintf('Transfer perigee speed        = %.6f km/s\n', v_perigee_trans)
fprintf('Transfer apogee speed         = %.6f km/s\n', v_apogee_trans)
fprintf('Approx Moon orbital speed     = %.6f km/s\n', v_moon_circ)
fprintf('Approx arrival v_infinity     = %.6f km/s\n', v_inf_moon)

%% ---------------- PART 6: LUNAR CAPTURE ESTIMATE ----------------
h_perilune = 100;               % km target perilune altitude
rp_moon    = RM + h_perilune;   % km

v_perilune_hyp = sqrt(v_inf_moon^2 + 2*muM/rp_moon);
v_circ_moon    = sqrt(muM/rp_moon);
delV_capture   = v_perilune_hyp - v_circ_moon;

fprintf('\n--- Lunar Capture Estimate ---\n')
fprintf('Target perilune altitude      = %.6f km\n', h_perilune)
fprintf('Hyperbolic perilune speed     = %.6f km/s\n', v_perilune_hyp)
fprintf('Circular lunar-orbit speed    = %.6f km/s\n', v_circ_moon)
fprintf('Capture Delta-V               = %.6f km/s\n', delV_capture)
fprintf('Approx total Delta-V          = %.6f km/s\n', delV1 + delV_capture)

%% ---------------- PART 7: PROPAGATE THE TLI TRANSFER ----------------
N_tli = 800;
tspan_tli = linspace(0, t_trans, N_tli);

r_tli = zeros(N_tli,3);
v_tli = zeros(N_tli,3);

for k = 1:N_tli
    [rtemp, vtemp] = UVars(rTLI_vec, vTLI_vec, tspan_tli(k));
    r_tli(k,:) = rtemp(:).';
    v_tli(k,:) = vtemp(:).';
end

r_tli_end = r_tli(end,:).';
v_tli_end = v_tli(end,:).';

fprintf('\nState at transfer apogee:\n')
fprintf('r_end = [%.6f %.6f %.6f] km\n', r_tli_end(1), r_tli_end(2), r_tli_end(3))
fprintf('v_end = [%.6f %.6f %.6f] km/s\n', v_tli_end(1), v_tli_end(2), v_tli_end(3))
fprintf('Final propagated radius       = %.6f km\n', norm(r_tli_end))
fprintf('Expected Moon distance        = %.6f km\n', r_moon)
fprintf('Radius difference             = %.6e km\n', norm(r_tli_end) - r_moon)

%% ---------------- PART 8: PLOTS ----------------
figure
plot3(r_tli(:,1), r_tli(:,2), r_tli(:,3), 'b-', 'LineWidth', 1.5)
hold on
plot3(rTLI_vec(1), rTLI_vec(2), rTLI_vec(3), 'go', 'MarkerFaceColor', 'g')
plot3(r_tli_end(1), r_tli_end(2), r_tli_end(3), 'mo', 'MarkerFaceColor', 'm')

[xE,yE,zE] = sphere(100);
surf(RE*xE, RE*yE, RE*zE, ...
    'FaceColor', [0.6 0.8 1.0], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.35)

grid on
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
legend('TLI transfer arc', 'TLI burn point', 'Transfer apogee', 'Earth', 'Location', 'best')
title('Hohmann-like Translunar Injection')
view(3)

figure
plot(tspan_tli/3600, vecnorm(r_tli,2,2), 'b-', 'LineWidth', 1.5)
grid on
xlabel('Time (hr)')
ylabel('|r| (km)')
        title('Radius Magnitude During Translunar Coast')

%% ---------------- PART 9: LAUNCH-TIME PHASING SWEEP ----------------
% Goal:
% Keep the same Hohmann-like transfer shape, but vary launch time to see
% when the spacecraft arrives closest to the Moon's actual position.

launch_date0 = datetime(2026,4,1,12,0,0);   % nominal launch date

% Sweep launch time around the nominal date
% Example: +/- 5 days in 1-hour increments
dt_hours = -5*24 : 1 : 5*24;
nSweep   = length(dt_hours);

miss_distance = zeros(nSweep,1);
launch_dates  = launch_date0 + hours(dt_hours);

r_moon_arr_all = zeros(nSweep,3);
r_sc_arr_all   = zeros(nSweep,3);

for j = 1:nSweep
    
    % Candidate launch and arrival times
    launch_date_j = launch_dates(j);
    arrival_date_j = launch_date_j + seconds(t_trans);
    jd_arrival_j = juliandate(arrival_date_j);
    
    % Moon state at arrival
    [r_moon_arr_j, v_moon_arr_j] = planetEphemeris(jd_arrival_j,'Earth','Moon');
    r_moon_arr_j = r_moon_arr_j(:);
    v_moon_arr_j = v_moon_arr_j(:); %#ok<NASGU>
    
    % For this Hohmann-style model, the departure state is still the same
    % parking-orbit state with the same tangential TLI burn
    r_depart_j = rTLI_vec;
    v_depart_j = vTLI_vec;
    
    % Propagate spacecraft to arrival using your existing UVars function
    [r_sc_arr_j, v_sc_arr_j] = UVars(r_depart_j, v_depart_j, t_trans); %#ok<ASGLU>
    r_sc_arr_j = r_sc_arr_j(:);
    
    % Miss distance
    miss_distance(j) = norm(r_sc_arr_j - r_moon_arr_j);
    
    % Store for plotting / reporting
    r_moon_arr_all(j,:) = r_moon_arr_j.';
    r_sc_arr_all(j,:)   = r_sc_arr_j.';
end

% Best launch time in the sweep
[miss_best, idx_best] = min(miss_distance);
best_launch_date  = launch_dates(idx_best);
best_arrival_date = best_launch_date + seconds(t_trans);

fprintf('\n--- Launch-Time Phasing Sweep ---\n')
fprintf('Nominal launch date           = %s\n', datestr(launch_date0))
fprintf('Best launch date in sweep     = %s\n', datestr(best_launch_date))
fprintf('Best arrival date             = %s\n', datestr(best_arrival_date))
fprintf('Minimum miss distance         = %.6f km\n', miss_best)

%% ---------------- PART 10: REPORT BEST-CASE ARRIVAL GEOMETRY ----------------
r_sc_best   = r_sc_arr_all(idx_best,:).';
r_moon_best = r_moon_arr_all(idx_best,:).';

fprintf('\nBest-case arrival vectors:\n')
fprintf('Spacecraft arrival position   = [%.6f %.6f %.6f] km\n', ...
    r_sc_best(1), r_sc_best(2), r_sc_best(3))
fprintf('Moon arrival position         = [%.6f %.6f %.6f] km\n', ...
    r_moon_best(1), r_moon_best(2), r_moon_best(3))
fprintf('Best miss distance            = %.6f km\n', norm(r_sc_best - r_moon_best))

%% ---------------- PART 11: PLOTS FOR PHASING SWEEP ----------------
figure
plot(dt_hours, miss_distance, 'b-', 'LineWidth', 1.5)
grid on
xlabel('Launch time offset from nominal (hr)')
ylabel('Miss distance at arrival (km)')
title('Moon Miss Distance vs Launch-Time Offset')

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
