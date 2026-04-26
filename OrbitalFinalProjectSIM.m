clc; clear; close all

%% ================================================================
%  BAREBONES EARTH-MOON MISSION SCRIPT WITH LUNAR GATEWAY RPO
% ================================================================

%% Constants
muE = 398600.4418;      % km^3/s^2
RE  = 6378.1363;        % km
muM = 4902.800066;      % km^3/s^2
RM  = 1737.4;           % km

launch_date0 = datetime(2026,4,1,12,0,0);

time_initial_LEO = launch_date0;
time_plane_change = launch_date0;   % simplified model: plane change occurs immediately

%% Mission Inputs
h_park = 185;           % km
inc_initial = 28.5;     % deg
inc_moon = 28.6;        % deg, simplified Earth-Moon plane inclination
RAAN = 0;               % deg
w = 0;                  % deg
theta = 0;              % deg
epark = 0;

rpark = RE + h_park;
r_moon_mean = 384400;   % km

h_lunar_final = 100;    % km final circular lunar orbit altitude
h_landing = 0;          % km landing altitude above lunar surface

rp_lunar_target = RM + h_lunar_final;
r_landing = RM + h_landing;

capture_radius = 10000; % km from Moon center where capture burn is attempted

% Lunar Gateway / RPO inputs
gateway_phase_deg = 20;           % Gateway is 20 deg ahead at circularization
rpo_phase_revs_min = 1;           % minimum phasing revolutions to attempt
rpo_min_altitude = 20;            % km, minimum allowable RPO phasing periapsis altitude
post_rendezvous_loiter_hr = 1.5;  % hr, coast with Gateway before deorbit

%% ================================================================
%  1. INITIAL LEO STATE
% ================================================================
hmag = sqrt(muE*rpark);

[r0, v0] = Perifocal2GE_simple(hmag, inc_initial, RAAN, epark, w, theta, muE);

%% ================================================================
%  2. SIMPLIFIED INCLINATION CHANGE BURN
% ================================================================
[~, v_after_plane] = Perifocal2GE_simple(hmag, inc_moon, RAAN, epark, w, theta, muE);

dV_plane_vec = v_after_plane - v0;
dV_plane = norm(dV_plane_vec);

r_current = r0(:);
v_current = v_after_plane(:);

%% ================================================================
%  3. COMPUTE HOHMANN-LIKE TLI BURN
% ================================================================
[dV_TLI, ~, t_trans] = Hohmann_simple(rpark, rpark, r_moon_mean, r_moon_mean, muE);

%% ================================================================
%  4. FIND BEST TLI TIME OFFSET USING SIMPLE FORWARD SWEEP
% ================================================================
dt_hours = 0:0.1:8;
miss = zeros(size(dt_hours));

for k = 1:length(dt_hours)

    dt_sec = dt_hours(k)*3600;

    [r_leo_k, v_leo_k] = propagateTwoBody(r_current, v_current, dt_sec, muE);

    vhat = v_leo_k/norm(v_leo_k);
    r_tli_k = r_leo_k;
    v_tli_k = v_leo_k + dV_TLI*vhat;

    launch_date_k = launch_date0 + seconds(dt_sec);
    jd_tli_k = juliandate(launch_date_k);

    [~, ~, r_end, ~, rMoon_end, ~] = propagateEarthMoon( ...
        r_tli_k, v_tli_k, jd_tli_k, t_trans, muE, muM, 600);

    miss(k) = norm(r_end - rMoon_end);

end

[~, idx_best] = min(miss);
dt_best = dt_hours(idx_best)*3600;

tli_date = launch_date0 + seconds(dt_best);
jd_tli = juliandate(tli_date);

time_TLI = tli_date;
t_TLI_mission = seconds(time_TLI - launch_date0);

%% ================================================================
%  5. PROPAGATE LEO TO TLI POINT, THEN DO TLI BURN
% ================================================================
[r_tli, v_before_tli] = propagateTwoBody(r_current, v_current, dt_best, muE);

vhat_tli = v_before_tli/norm(v_before_tli);
v_after_tli = v_before_tli + dV_TLI*vhat_tli;

%% ================================================================
%  6. PROPAGATE TLI TO NEAR-MOON CAPTURE LOCATION
% ================================================================
[t_TLI, y_TLI, r_capture_ECI, v_capture_ECI, rMoon_capture, vMoon_capture, t_capture] = ...
    propagateToMoonCapture(r_tli, v_after_tli, jd_tli, t_trans, ...
                           capture_radius, muE, muM, 2000);

time_capture = time_TLI + seconds(t_capture);
t_capture_mission = seconds(time_capture - launch_date0);

%% ================================================================
%  7. DO LUNAR CAPTURE BURN
% ================================================================
r_rel_capture = r_capture_ECI - rMoon_capture;
v_rel_capture_before = v_capture_ECI - vMoon_capture;

r_cap = norm(r_rel_capture);

% Simplified captured ellipse:
% Treat the capture point as apolune and target a low-lunar-orbit perilune.
ra_lunar = r_cap;
rp_lunar = rp_lunar_target;
a_lunar = 0.5*(ra_lunar + rp_lunar);

% Speed required at apolune of this lunar ellipse
v_capture_target = sqrt(muM*(2/ra_lunar - 1/a_lunar));

% Force the post-capture velocity to be tangential at apolune.
rhat_cap = r_rel_capture / norm(r_rel_capture);

hhat_cap = cross(r_rel_capture, v_rel_capture_before);

if norm(hhat_cap) < 1e-12
    hhat_cap = cross(rhat_cap, [0;0;1]);
    if norm(hhat_cap) < 1e-12
        hhat_cap = cross(rhat_cap, [0;1;0]);
    end
end

hhat_cap = hhat_cap / norm(hhat_cap);

vhat_tan_cap = cross(hhat_cap, rhat_cap);
vhat_tan_cap = vhat_tan_cap / norm(vhat_tan_cap);

% Post-capture Moon-relative velocity
v_rel_capture_after = v_capture_target * vhat_tan_cap;

dV_capture_vec = v_rel_capture_after - v_rel_capture_before;
dV_capture = norm(dV_capture_vec);

%% ================================================================
%  8. PROPAGATE CAPTURED LUNAR ELLIPSE TO PERILUNE
% ================================================================
% Since the capture point is constructed as apolune, perilune occurs after
% half of the captured lunar ellipse period.

t_to_perilune = pi*sqrt(a_lunar^3/muM);

[t_lunar_ellipse, y_lunar_ellipse] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,t_to_perilune,800), ...
    [r_rel_capture; v_rel_capture_after], ...
    odeset('RelTol',1e-11,'AbsTol',1e-11));

time_perilune = time_capture + seconds(t_to_perilune);
time_circularization = time_perilune;

t_perilune_mission = seconds(time_perilune - launch_date0);
t_circularization_mission = t_perilune_mission;

r_peri = y_lunar_ellipse(end,1:3).';
v_peri_before = y_lunar_ellipse(end,4:6).';

r_peri_mag = norm(r_peri);

%% ================================================================
%  9. DO LUNAR CIRCULARIZATION BURN
% ================================================================
v_circ_mag = sqrt(muM/r_peri_mag);

% Make the circularized velocity exactly tangential at perilune.
rhat_peri = r_peri / norm(r_peri);
hhat_peri = cross(r_peri, v_peri_before);
hhat_peri = hhat_peri / norm(hhat_peri);

vhat_peri_tan = cross(hhat_peri, rhat_peri);
vhat_peri_tan = vhat_peri_tan / norm(vhat_peri_tan);

v_peri_after = v_circ_mag * vhat_peri_tan;

dV_circ_vec = v_peri_after - v_peri_before;
dV_circ = norm(dV_circ_vec);

%% ================================================================
%  10. PROPAGATE FINAL CIRCULAR LUNAR ORBIT REFERENCE
% ================================================================
T_circ = 2*pi*sqrt(r_peri_mag^3/muM);

[t_lunar_circ, y_lunar_circ] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,T_circ,800), ...
    [r_peri; v_peri_after], ...
    odeset('RelTol',1e-12,'AbsTol',1e-12));

time_final_orbit_complete = time_circularization + seconds(T_circ);
t_final_orbit_complete_mission = seconds(time_final_orbit_complete - launch_date0);

%% ================================================================
%  10A. RPO PHASING RENDEZVOUS WITH LUNAR GATEWAY
% ================================================================
% Gateway is in the same 100 km circular lunar orbit and is initially
% 20 deg ahead of the spacecraft at circularization.
%
% Since the Gateway is ahead, the spacecraft drops to a slightly lower,
% faster phasing orbit. The number of phasing revolutions is increased
% automatically until the phasing periapsis stays above the requested
% minimum altitude.

r_rpo_start = r_peri;
v_rpo_circ_start = v_peri_after;

r_circ_rpo = norm(r_rpo_start);
v_circ_rpo = norm(v_rpo_circ_start);

hhat_orbit = cross(r_rpo_start, v_rpo_circ_start);
hhat_orbit = hhat_orbit / norm(hhat_orbit);

% Gateway initial state, 20 deg ahead in the direction of orbital motion
r_gateway0 = rotateAboutAxis(r_rpo_start, hhat_orbit, deg2rad(gateway_phase_deg));
v_gateway0 = rotateAboutAxis(v_rpo_circ_start, hhat_orbit, deg2rad(gateway_phase_deg));

% Circular orbit period and mean motion
T_gateway = 2*pi*sqrt(r_circ_rpo^3/muM);

% Gateway is ahead by phi. Spacecraft must gain phi over N phasing revs.
phi = deg2rad(gateway_phase_deg);

% Choose valid number of phasing revolutions
rpo_phase_revs = rpo_phase_revs_min;

while true

    T_phase = T_gateway * (1 - phi/(2*pi*rpo_phase_revs));
    a_phase = (muM*(T_phase/(2*pi))^2)^(1/3);

    ra_phase = r_circ_rpo;
    rp_phase = 2*a_phase - ra_phase;

    if rp_phase > RM + rpo_min_altitude
        break
    end

    rpo_phase_revs = rpo_phase_revs + 1;

    if rpo_phase_revs > 50
        error('No valid RPO phasing orbit found. Increase final lunar orbit altitude or lower rpo_min_altitude.')
    end

end

% Speed at apolune of the lower phasing orbit
v_phase_apo = sqrt(muM*(2/ra_phase - 1/a_phase));

% First RPO burn: reduce speed into lower phasing orbit
vhat_rpo = v_rpo_circ_start / norm(v_rpo_circ_start);
v_rpo_depart = v_phase_apo * vhat_rpo;

dV_rpo_depart_vec = v_rpo_depart - v_rpo_circ_start;
dV_rpo_depart = norm(dV_rpo_depart_vec);

% Total RPO phasing time
t_rpo = rpo_phase_revs * T_phase;

time_rpo_start = time_circularization;
t_rpo_start_mission = seconds(time_rpo_start - launch_date0);

% Propagate spacecraft on phasing orbit
[t_rpo_phase, y_rpo_phase] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,t_rpo,1200), ...
    [r_rpo_start; v_rpo_depart], ...
    odeset('RelTol',1e-12,'AbsTol',1e-12));

% Propagate Gateway on circular orbit
[t_gateway, y_gateway] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,t_rpo,1200), ...
    [r_gateway0; v_gateway0], ...
    odeset('RelTol',1e-12,'AbsTol',1e-12));

% End states
r_rpo_numerical_end = y_rpo_phase(end,1:3).';
v_rpo_numerical_end = y_rpo_phase(end,4:6).';

r_gateway_rendezvous = y_gateway(end,1:3).';
v_gateway_rendezvous = y_gateway(end,4:6).';

% Use Gateway endpoint as the rendezvous point
r_rendezvous = r_gateway_rendezvous;

% Second RPO burn: circularize/match Gateway velocity
v_rendezvous_after = v_gateway_rendezvous;

dV_rpo_arrive_vec = v_rendezvous_after - v_rpo_numerical_end;
dV_rpo_arrive = norm(dV_rpo_arrive_vec);

dV_RPO_total = dV_rpo_depart + dV_rpo_arrive;

time_rendezvous = time_rpo_start + seconds(t_rpo);
t_rendezvous_mission = seconds(time_rendezvous - launch_date0);

% Diagnostic errors
rpo_position_error = norm(r_rpo_numerical_end - r_gateway_rendezvous);
rpo_radius_error = norm(r_rendezvous) - (RM + h_lunar_final);

%% ================================================================
%  10B. POST-RENDEZVOUS LOITER WITH GATEWAY
% ================================================================
% After rendezvous, spacecraft has matched Gateway velocity and coasts with
% Gateway in the 100 km circular lunar orbit before deorbit.

t_loiter = post_rendezvous_loiter_hr * 3600;

[t_loiter_hist, y_loiter] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,t_loiter,600), ...
    [r_rendezvous; v_rendezvous_after], ...
    odeset('RelTol',1e-12,'AbsTol',1e-12));

r_after_loiter = y_loiter(end,1:3).';
v_after_loiter = y_loiter(end,4:6).';

time_deorbit = time_rendezvous + seconds(t_loiter);
t_deorbit_mission = seconds(time_deorbit - launch_date0);

%% ================================================================
%  10C. DEORBIT AND LANDING BURN AFTER GATEWAY LOITER
% ================================================================
% Start deorbit after loitering with Gateway.

r_deorbit = r_after_loiter;
v_circ_before_deorbit = v_after_loiter;

r_deorbit_mag = norm(r_deorbit);

ra_descent = r_deorbit_mag;
rp_descent = r_landing;
a_descent = 0.5*(ra_descent + rp_descent);

% Speed at apolune of the descent ellipse
v_deorbit_mag = sqrt(muM*(2/ra_descent - 1/a_descent));

% Keep velocity tangential and reduce speed
vhat_deorbit = v_circ_before_deorbit / norm(v_circ_before_deorbit);
v_after_deorbit = v_deorbit_mag * vhat_deorbit;

dV_deorbit_vec = v_after_deorbit - v_circ_before_deorbit;
dV_deorbit = norm(dV_deorbit_vec);

% Coast from apolune of the descent ellipse to the lunar surface
t_to_landing = pi*sqrt(a_descent^3/muM);

[t_landing, y_landing] = ode45( ...
    @(t,y) stateTwoBody(t,y,muM), ...
    linspace(0,t_to_landing,800), ...
    [r_deorbit; v_after_deorbit], ...
    odeset('RelTol',1e-12,'AbsTol',1e-12));

time_landing = time_deorbit + seconds(t_to_landing);
t_landing_mission = seconds(time_landing - launch_date0);

r_land = y_landing(end,1:3).';
v_land_before = y_landing(end,4:6).';

% Landing burn: cancel all Moon-relative velocity at the surface
v_land_after = [0;0;0];

dV_landing_vec = v_land_after - v_land_before;
dV_landing = norm(dV_landing_vec);

%% ================================================================
%  11. PLOTS
% ================================================================

%% Parking orbit timing
T_park = 2*pi*sqrt(rpark^3/muE);
num_parking_orbits = dt_best/T_park;

% Original pre-plane-change parking orbit for reference
[~, y_initial_orbit] = ode45( ...
    @(t,y) stateTwoBody(t,y,muE), ...
    linspace(0,T_park,700), ...
    [r0(:); v0(:)], ...
    odeset('RelTol',1e-11,'AbsTol',1e-11));

% Actual post-plane-change coast before TLI
if dt_best > 1e-9
    [~, y_pc_coast] = ode45( ...
        @(t,y) stateTwoBody(t,y,muE), ...
        linspace(0,dt_best,1000), ...
        [r_current(:); v_current(:)], ...
        odeset('RelTol',1e-11,'AbsTol',1e-11));
else
    y_pc_coast = [r_current(:); v_current(:)].';
end

% One full post-plane-change parking orbit for visual reference
[~, y_post_plane_ref] = ode45( ...
    @(t,y) stateTwoBody(t,y,muE), ...
    linspace(0,T_park,700), ...
    [r_current(:); v_current(:)], ...
    odeset('RelTol',1e-11,'AbsTol',1e-11));

%% 1. Earth-centered parking orbit / plane-change / TLI plot
figure
hold on
grid on
axis equal
xlabel('x_E (km)')
ylabel('y_E (km)')
zlabel('z_E (km)')
title(sprintf('LEO Parking Orbit Before TLI: %.3f Parking Orbits Before TLI', ...
      num_parking_orbits))
view(3)

[xE,yE,zE] = sphere(60);
surf(RE*xE, RE*yE, RE*zE, ...
    'FaceAlpha',0.25, ...
    'EdgeColor','none')

plot3(y_initial_orbit(:,1), y_initial_orbit(:,2), y_initial_orbit(:,3), ...
      'k--', 'LineWidth', 1.0)

plot3(y_post_plane_ref(:,1), y_post_plane_ref(:,2), y_post_plane_ref(:,3), ...
      'b:', 'LineWidth', 1.0)

plot3(y_pc_coast(:,1), y_pc_coast(:,2), y_pc_coast(:,3), ...
      'b-', 'LineWidth', 1.6)

plot3(r_current(1), r_current(2), r_current(3), ...
      'mo', 'MarkerFaceColor','m', 'MarkerSize', 7)

plot3(r_tli(1), r_tli(2), r_tli(3), ...
      'ro', 'MarkerFaceColor','r', 'MarkerSize', 7)

legend('Earth', ...
       'Original LEO before plane change', ...
       'Post-plane-change LEO reference', ...
       'Actual coast before TLI', ...
       'Plane change burn', ...
       'TLI burn', ...
       'Location','best')

%% 2. Full Earth-centered transfer to Moon
figure
hold on
grid on
axis equal
xlabel('x_E (km)')
ylabel('y_E (km)')
zlabel('z_E (km)')
title('Earth-Centered Translunar Trajectory')
view(3)

surf(RE*xE, RE*yE, RE*zE, ...
    'FaceAlpha',0.25, ...
    'EdgeColor','none')

plot3(y_pc_coast(:,1), y_pc_coast(:,2), y_pc_coast(:,3), ...
      'b-', 'LineWidth', 1.3)

plot3(y_TLI(:,1), y_TLI(:,2), y_TLI(:,3), ...
      'r-', 'LineWidth', 1.3)

plot3(rMoon_capture(1), rMoon_capture(2), rMoon_capture(3), ...
      'ko', 'MarkerFaceColor','k')

plot3(r_current(1), r_current(2), r_current(3), ...
      'mo', 'MarkerFaceColor','m')

plot3(r_tli(1), r_tli(2), r_tli(3), ...
      'ro', 'MarkerFaceColor','r')

legend('Earth', ...
       'Parking orbit coast', ...
       'TLI trajectory', ...
       'Moon at capture', ...
       'Plane change burn', ...
       'TLI burn', ...
       'Location','best')

%% 3. Moon-centered capture, circularization, RPO, loiter, and landing plot
figure
hold on
grid on
axis equal
xlabel('x_M (km)')
ylabel('y_M (km)')
zlabel('z_M (km)')
title('Moon-Centered Capture, RPO, Loiter, and Landing')
view(3)

[xM,yM,zM] = sphere(60);
surf(RM*xM, RM*yM, RM*zM, ...
    'FaceAlpha',0.35, ...
    'EdgeColor','none')

plot3(y_lunar_ellipse(:,1), y_lunar_ellipse(:,2), y_lunar_ellipse(:,3), ...
      'm-', 'LineWidth', 1.3)

plot3(y_lunar_circ(:,1), y_lunar_circ(:,2), y_lunar_circ(:,3), ...
      'g-', 'LineWidth', 1.0)

plot3(y_rpo_phase(:,1), y_rpo_phase(:,2), y_rpo_phase(:,3), ...
      'b-', 'LineWidth', 1.6)

plot3(y_gateway(:,1), y_gateway(:,2), y_gateway(:,3), ...
      'k--', 'LineWidth', 1.1)

plot3(y_loiter(:,1), y_loiter(:,2), y_loiter(:,3), ...
      'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 1.4)

plot3(y_landing(:,1), y_landing(:,2), y_landing(:,3), ...
      'c-', 'LineWidth', 1.5)

plot3(r_rel_capture(1), r_rel_capture(2), r_rel_capture(3), ...
      'bo', 'MarkerFaceColor','b')

plot3(r_peri(1), r_peri(2), r_peri(3), ...
      'ro', 'MarkerFaceColor','r')

plot3(r_gateway0(1), r_gateway0(2), r_gateway0(3), ...
      'kd', 'MarkerFaceColor','k', 'MarkerSize', 7)

plot3(r_rendezvous(1), r_rendezvous(2), r_rendezvous(3), ...
      'gs', 'MarkerFaceColor','g', 'MarkerSize', 7)

plot3(r_deorbit(1), r_deorbit(2), r_deorbit(3), ...
      'ko', 'MarkerFaceColor','k', 'MarkerSize', 7)

plot3(r_land(1), r_land(2), r_land(3), ...
      'ys', 'MarkerFaceColor','y', 'MarkerSize', 7)

moon_view_radius = max([capture_radius, ...
                        1.2*max(vecnorm(y_lunar_ellipse(:,1:3),2,2)), ...
                        1.2*max(vecnorm(y_rpo_phase(:,1:3),2,2)), ...
                        1.2*max(vecnorm(y_loiter(:,1:3),2,2)), ...
                        1.2*max(vecnorm(y_landing(:,1:3),2,2))]);

xlim([-moon_view_radius moon_view_radius])
ylim([-moon_view_radius moon_view_radius])
zlim([-moon_view_radius moon_view_radius])

legend('Moon actual size', ...
       'Captured lunar ellipse', ...
       '100 km circular orbit reference', ...
       'RPO phasing trajectory', ...
       'Gateway trajectory', ...
       'Post-rendezvous loiter', ...
       'Descent trajectory', ...
       'Capture burn', ...
       'Circularization / RPO start', ...
       'Gateway initial position', ...
       'Rendezvous', ...
       'Deorbit burn after loiter', ...
       'Landing burn / touchdown', ...
       'Location','best')

%% 4. Orbital-plane RPO rendezvous plot
figure
hold on
grid on
axis equal
xlabel('In-plane x_M (km)')
ylabel('In-plane y_M (km)')
title('Lunar Gateway RPO Rendezvous: Orbital Plane View')

% Build orbital-plane basis
e1 = r_rpo_start / norm(r_rpo_start);
e2 = cross(hhat_orbit, e1);
e2 = e2 / norm(e2);

% Helper projection
project2D = @(r) [dot(r(:),e1), dot(r(:),e2)];

% Moon and 100 km circular orbit in the orbital plane
ang = linspace(0,2*pi,500);

moon_x = RM*cos(ang);
moon_y = RM*sin(ang);

orbit_radius = RM + h_lunar_final;
orbit_x = orbit_radius*cos(ang);
orbit_y = orbit_radius*sin(ang);

plot(moon_x, moon_y, 'k-', 'LineWidth', 1.2)
plot(orbit_x, orbit_y, 'g:', 'LineWidth', 1.5)

% Project spacecraft RPO trajectory
rpo_xy = zeros(length(t_rpo_phase),2);
for k = 1:length(t_rpo_phase)
    rpo_xy(k,:) = project2D(y_rpo_phase(k,1:3).');
end

% Project Gateway trajectory
gateway_xy = zeros(length(t_gateway),2);
for k = 1:length(t_gateway)
    gateway_xy(k,:) = project2D(y_gateway(k,1:3).');
end

% Project loiter trajectory
loiter_xy = zeros(length(t_loiter_hist),2);
for k = 1:length(t_loiter_hist)
    loiter_xy(k,:) = project2D(y_loiter(k,1:3).');
end

plot(rpo_xy(:,1), rpo_xy(:,2), 'b-', 'LineWidth', 1.6)
plot(gateway_xy(:,1), gateway_xy(:,2), 'k--', 'LineWidth', 1.2)
plot(loiter_xy(:,1), loiter_xy(:,2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.4)

r_rpo_start_xy = project2D(r_rpo_start);
r_gateway0_xy = project2D(r_gateway0);
r_rendezvous_xy = project2D(r_rendezvous);
r_rpo_end_xy = project2D(r_rpo_numerical_end);
r_deorbit_xy = project2D(r_deorbit);

plot(r_rpo_start_xy(1), r_rpo_start_xy(2), ...
     'ro', 'MarkerFaceColor','r', 'MarkerSize', 7)

plot(r_gateway0_xy(1), r_gateway0_xy(2), ...
     'kd', 'MarkerFaceColor','k', 'MarkerSize', 7)

plot(r_rendezvous_xy(1), r_rendezvous_xy(2), ...
     'gs', 'MarkerFaceColor','g', 'MarkerSize', 8)

plot(r_rpo_end_xy(1), r_rpo_end_xy(2), ...
     'bo', 'MarkerFaceColor','b', 'MarkerSize', 6)

plot(r_deorbit_xy(1), r_deorbit_xy(2), ...
     'ko', 'MarkerFaceColor','k', 'MarkerSize', 7)

legend('Moon radius', ...
       '100 km lunar orbit', ...
       'Spacecraft RPO phasing trajectory', ...
       'Gateway trajectory', ...
       'Post-rendezvous loiter', ...
       'RPO departure', ...
       'Gateway initially 20 deg ahead', ...
       'Rendezvous point on Gateway orbit', ...
       'Propagated spacecraft endpoint', ...
       'Deorbit after loiter', ...
       'Location','best')

view(2)

%% 5. Zoomed Moon-centered landing burn propagation plot
figure
hold on
grid on
axis equal
xlabel('x_M (km)')
ylabel('y_M (km)')
zlabel('z_M (km)')
title('Zoomed Lunar Deorbit and Landing Burn Propagation')
view(3)

[xM2,yM2,zM2] = sphere(80);
surf(RM*xM2, RM*yM2, RM*zM2, ...
    'FaceAlpha',0.35, ...
    'EdgeColor','none')

plot3(y_landing(:,1), y_landing(:,2), y_landing(:,3), ...
      'c-', 'LineWidth', 1.8)

plot3(r_deorbit(1), r_deorbit(2), r_deorbit(3), ...
      'ko', 'MarkerFaceColor','k', 'MarkerSize', 7)

plot3(r_land(1), r_land(2), r_land(3), ...
      'ys', 'MarkerFaceColor','y', 'MarkerSize', 8)

plot3([0 r_land(1)], [0 r_land(2)], [0 r_land(3)], ...
      'k--', 'LineWidth', 0.8)

landing_zoom = 500; % km around touchdown point

xlim([r_land(1)-landing_zoom, r_land(1)+landing_zoom])
ylim([r_land(2)-landing_zoom, r_land(2)+landing_zoom])
zlim([r_land(3)-landing_zoom, r_land(3)+landing_zoom])

legend('Moon actual size', ...
       'Descent trajectory', ...
       'Deorbit burn after loiter', ...
       'Landing burn / touchdown', ...
       'Landing radial line', ...
       'Location','best')

%% ================================================================
%  12. SIMPLE NUMERIC SUMMARY
% ================================================================
dV_total = dV_plane + dV_TLI + dV_capture + dV_circ + ...
           dV_rpo_depart + dV_rpo_arrive + dV_deorbit + dV_landing;

fprintf('\nDelta-V Summary:\n')
fprintf('Plane change          = %.6f km/s\n', dV_plane)
fprintf('TLI                   = %.6f km/s\n', dV_TLI)
fprintf('Lunar capture         = %.6f km/s\n', dV_capture)
fprintf('Circularization       = %.6f km/s\n', dV_circ)
fprintf('RPO departure         = %.6f km/s\n', dV_rpo_depart)
fprintf('RPO rendezvous        = %.6f km/s\n', dV_rpo_arrive)
fprintf('Deorbit               = %.6f km/s\n', dV_deorbit)
fprintf('Landing burn          = %.6f km/s\n', dV_landing)
fprintf('Total                 = %.6f km/s\n', dV_total)

fprintf('\nMission Timeline:\n')
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Initial LEO / launch', 0, datestr(time_initial_LEO))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Plane change burn', 0, datestr(time_plane_change))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'TLI burn', t_TLI_mission/3600, datestr(time_TLI))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Lunar capture burn', t_capture_mission/3600, datestr(time_capture))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Perilune reached', t_perilune_mission/3600, datestr(time_perilune))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Circularization burn', t_circularization_mission/3600, datestr(time_circularization))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'RPO departure burn', t_rpo_start_mission/3600, datestr(time_rpo_start))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Gateway rendezvous', t_rendezvous_mission/3600, datestr(time_rendezvous))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'End Gateway loiter', t_deorbit_mission/3600, datestr(time_deorbit))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Deorbit burn', t_deorbit_mission/3600, datestr(time_deorbit))
fprintf('%-30s  MET = %10.3f hr   Date = %s\n', ...
    'Landing burn / touchdown', t_landing_mission/3600, datestr(time_landing))

fprintf('\nRPO Check:\n')
fprintf('Gateway initial lead angle             = %.6f deg\n', gateway_phase_deg)
fprintf('RPO phasing revolutions                = %d\n', rpo_phase_revs)
fprintf('RPO phasing time                       = %.6f hr\n', t_rpo/3600)
fprintf('RPO phasing periapsis altitude         = %.6f km\n', rp_phase - RM)
fprintf('Post-rendezvous loiter time            = %.6f hr\n', post_rendezvous_loiter_hr)
fprintf('Rendezvous radius                      = %.6f km\n', norm(r_rendezvous))
fprintf('Rendezvous altitude                    = %.6f km\n', norm(r_rendezvous) - RM)
fprintf('Rendezvous radius error from 100 km orbit = %.6e km\n', rpo_radius_error)
fprintf('Spacecraft endpoint error              = %.6f km\n', rpo_position_error)
fprintf('Rendezvous velocity before arrival burn= %.6f km/s\n', norm(v_rpo_numerical_end - v_gateway_rendezvous))

fprintf('\nLunar Orbit and Landing Check:\n')
fprintf('Target perilune altitude               = %.6f km\n', h_lunar_final)
fprintf('Actual perilune altitude               = %.6f km\n', r_peri_mag - RM)
fprintf('Final circular altitude                = %.6f km\n', mean(vecnorm(y_lunar_circ(:,1:3),2,2)) - RM)
fprintf('Landing altitude                       = %.6f km\n', norm(r_land) - RM)
fprintf('Landing speed before burn              = %.6f km/s\n', norm(v_land_before))

%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function [r_vec, v_vec] = Perifocal2GE_simple(h, inc_deg, RAAN_deg, e, w_deg, theta_deg, mu)

    i = deg2rad(inc_deg);
    RAAN = deg2rad(RAAN_deg);
    w = deg2rad(w_deg);
    theta = deg2rad(theta_deg);

    r_pf = (h^2/mu)/(1 + e*cos(theta))*[cos(theta); sin(theta); 0];
    v_pf = (mu/h)*[-sin(theta); e + cos(theta); 0];

    R3_W = [ cos(RAAN) -sin(RAAN) 0;
             sin(RAAN)  cos(RAAN) 0;
             0          0         1];

    R1_i = [1 0       0;
            0 cos(i) -sin(i);
            0 sin(i)  cos(i)];

    R3_w = [ cos(w) -sin(w) 0;
             sin(w)  cos(w) 0;
             0       0      1];

    Q = R3_W*R1_i*R3_w;

    r_vec = Q*r_pf;
    v_vec = Q*v_pf;

end

function [dV1, dV2, t_trans] = Hohmann_simple(rA, rAp, rB, rBp, mu)

    h1 = sqrt(2*mu)*sqrt(rA*rAp/(rA+rAp));
    h2 = sqrt(2*mu)*sqrt(rB*rBp/(rB+rBp));
    h3 = sqrt(2*mu)*sqrt(rA*rB/(rA+rB));

    vA1 = h1/rA;
    vA3 = h3/rA;
    dV1 = vA3 - vA1;

    vB3 = h3/rB;
    vB2 = h2/rB;
    dV2 = vB2 - vB3;

    a_trans = 0.5*(rA+rB);
    t_trans = pi*sqrt(a_trans^3/mu);

end

function [rout, vout] = propagateTwoBody(r0, v0, tf, mu)

    r0 = r0(:);
    v0 = v0(:);

    if abs(tf) < 1e-12
        rout = r0;
        vout = v0;
        return
    end

    [~, y] = ode45(@(t,y) stateTwoBody(t,y,mu), ...
                   [0 tf], ...
                   [r0; v0], ...
                   odeset('RelTol',1e-11,'AbsTol',1e-11));

    rout = y(end,1:3).';
    vout = y(end,4:6).';

end

function [t_hist, y_hist, r_end, v_end, rMoon_end, vMoon_end] = ...
    propagateEarthMoon(r0, v0, jd0, tf, muE, muM, N)

    [t_hist, y_hist] = ode45(@(t,y) stateEarthMoon(t,y,muE,muM,jd0), ...
                             linspace(0,tf,N), ...
                             [r0(:); v0(:)], ...
                             odeset('RelTol',1e-10,'AbsTol',1e-10));

    r_end = y_hist(end,1:3).';
    v_end = y_hist(end,4:6).';

    jd_end = jd0 + t_hist(end)/86400;
    [rMoon_end, vMoon_end] = planetEphemeris(jd_end,'Earth','Moon');

    rMoon_end = rMoon_end(:);
    vMoon_end = vMoon_end(:);

end

function [t_hist, y_hist, r_cap, v_cap, rMoon_cap, vMoon_cap, t_cap] = ...
    propagateToMoonCapture(r0, v0, jd0, tf, capture_radius, muE, muM, N)

    opts = odeset('RelTol',1e-10, ...
                  'AbsTol',1e-10, ...
                  'Events', @(t,y) moonCaptureEvent(t,y,jd0,capture_radius));

    [t_hist, y_hist, t_event, y_event] = ode45( ...
        @(t,y) stateEarthMoon(t,y,muE,muM,jd0), ...
        linspace(0,tf,N), ...
        [r0(:); v0(:)], ...
        opts);

    if isempty(t_event)

        r_all = y_hist(:,1:3);
        sep = zeros(size(t_hist));

        for k = 1:length(t_hist)
            jd = jd0 + t_hist(k)/86400;
            rMoon = planetEphemeris(jd,'Earth','Moon');
            sep(k) = norm(r_all(k,:).' - rMoon(:));
        end

        [~, idx] = min(sep);
        t_cap = t_hist(idx);
        y_cap = y_hist(idx,:).';

    else

        t_cap = t_event(end);
        y_cap = y_event(end,:).';

    end

    r_cap = y_cap(1:3);
    v_cap = y_cap(4:6);

    jd_cap = jd0 + t_cap/86400;
    [rMoon_cap, vMoon_cap] = planetEphemeris(jd_cap,'Earth','Moon');

    rMoon_cap = rMoon_cap(:);
    vMoon_cap = vMoon_cap(:);

end

function dydt = stateTwoBody(~, y, mu)

    r = y(1:3);
    v = y(4:6);

    a = -mu*r/norm(r)^3;

    dydt = [v; a];

end

function dydt = stateEarthMoon(t, y, muE, muM, jd0)

    r = y(1:3);
    v = y(4:6);

    jd = jd0 + t/86400;

    [rMoon, ~] = planetEphemeris(jd,'Earth','Moon');
    rMoon = rMoon(:);

    aE = -muE*r/norm(r)^3;
    aM = -muM*(r - rMoon)/norm(r - rMoon)^3;

    dydt = [v; aE + aM];

end

function [value, isterminal, direction] = moonCaptureEvent(t, y, jd0, capture_radius)

    r = y(1:3);

    jd = jd0 + t/86400;
    rMoon = planetEphemeris(jd,'Earth','Moon');
    rMoon = rMoon(:);

    value = norm(r - rMoon) - capture_radius;

    isterminal = 1;
    direction = -1;

end

function r_rot = rotateAboutAxis(r_vec, axis_vec, angle_rad)

    r_vec = r_vec(:);
    axis_vec = axis_vec(:)/norm(axis_vec);

    r_rot = r_vec*cos(angle_rad) + ...
            cross(axis_vec, r_vec)*sin(angle_rad) + ...
            axis_vec*dot(axis_vec, r_vec)*(1 - cos(angle_rad));

end
