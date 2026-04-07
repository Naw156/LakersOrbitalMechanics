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

#### NATHANS PART TO GET TO MOON!!!!###################################
clc
clear all

% 1. Set your Launch Date (Straight into the code)
launch_date = datetime(2026, 4, 1, 12, 0, 0); % April 7, 2026, at Noon

% 2. Define your travel time in seconds (e.g., 5 days)
t = 432000; 

arrival_jd = juliandate(launch_date);

% 4. Get the Moon's Position (rB_vec) from the Ephemeris
[rB_vec, v_moon] = planetEphemeris(arrival_jd, 'Earth', 'Moon');

moonr = rB_vec;
moonv = v_moon;
mu = 4.035*10^14;
ri = [6563.136, 0, 0];
vi = [0,6.849, 3.719];
[rB_vec,vA_vec,vB_vec,delVA] = Chase(ri, vi, moonr, moonv,t, mu, 'p');
