%% ========================================================================
%  EEEN415 T5 Presentation — Inverted Pendulum on a Cart
%  Section 1: Plant model + open-loop analysis
% =========================================================================

%% Plant parameters
M = 0.5;    % mass of cart            [kg]
m = 0.2;    % mass of pendulum        [kg]
b = 0.1;    % cart friction coeff     [N/m/s]
I = 0.006;  % pendulum inertia        [kg.m^2]
g = 9.8;    % gravity                 [m/s^2]
l = 0.3;    % length to pend. CoM     [m]

p = I*(M+m)+M*m*l^2;   % common denominator for A,B

%% State-space model:  states = [x, x_dot, phi, phi_dot]
A = [0      1              0           0;
    0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
    0      0              0           1;
    0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
    (I+m*l^2)/p;
    0;
    m*l/p];
C = [1 0 0 0;
    0 0 1 0];
D = [0;
    0];

states  = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs  = {'u'};
outputs = {'x'; 'phi'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs)

%% ------------------------------------------------------------------------
%  Controllability and open-loop poles
% -------------------------------------------------------------------------
Co = ctrb(A,B);
r  = rank(Co);
fprintf('Controllability matrix rank = %d of %d states\n', r, size(A,1));
if r == size(A,1)
    fprintf('System is controllable -> arbitrary pole placement is possible.\n');
else
    fprintf('System is NOT fully controllable.\n');
end

p_ol = eig(A);
disp('Open-loop poles:');
disp(p_ol);

%% Pole plot on the s-plane
figure; hold on; grid on;
plot(real(p_ol), imag(p_ol), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1);   % imaginary axis = stability boundary
yline(0, 'k-');
xlabel('Real axis  (\sigma)');
ylabel('Imaginary axis  (j\omega)');
title('Open-loop poles of the inverted pendulum');
legend('Open-loop poles', 'Location', 'best');
xlim([min(real(p_ol))-1, max(real(p_ol))+1]);
axis equal;

%% ------------------------------------------------------------------------
%  Open-loop regulator test: release from an initial tilt, no control
%  (shows the uncontrolled system diverges -- driven by the RHP pole)
% -------------------------------------------------------------------------
t  = 0:0.001:5;                    % fine step so the validity crossing is accurate
x0 = [0; 0; 0.1; 0];               % 0.1 rad initial tilt, cart at rest at origin

[y,t] = initial(sys_ss, x0, t);    % u = 0 throughout

phi_limit = 0.35;                  % 20 deg -- small-angle / spec validity limit
i_lim     = find(y(:,2) > phi_limit, 1, 'first');
t_lim     = t(i_lim);
fprintf('Uncontrolled: pendulum exceeds %.2f rad (20 deg) at t = %.3f s (x = %.4f m)\n', ...
        phi_limit, t_lim, y(i_lim,1));

figure;

% --- Cart position ---
subplot(2,1,1); hold on; grid on;
plot(t, y(:,1), 'LineWidth', 1.4);
xline(t_lim, 'k--', 'LineWidth', 1.2);
ylabel('Cart position x (m)');
title('Uncontrolled response to a 0.1 rad initial tilt');
xlim([0 0.5]); ylim([0 0.07]);
legend('x', 'Linearisation validity limit', 'Location', 'northwest');

% --- Pendulum angle ---
subplot(2,1,2); hold on; grid on;
plot(t, y(:,2), 'r', 'LineWidth', 1.4);
yline(phi_limit, 'k:',  'LineWidth', 1.2);
xline(t_lim,     'k--', 'LineWidth', 1.2);
ylabel('Pendulum angle \phi (rad)');
xlabel('Time (s)');
xlim([0 0.5]); ylim([0 0.9]);
legend('\phi', '20\circ limit (0.35 rad)', ...
       sprintf('Limit crossed at t = %.2f s', t_lim), 'Location', 'northwest');


%% ========================================================================
%  Section 2: Pole placement design
% =========================================================================

%% Desired pole locations from the design specifications
zeta = 0.7;               % damping ratio  -> low overshoot
wn   = 5;                 % natural freq   -> fast enough to catch the fall
%                   inside the 0.35 s validity window

s_dom = [-zeta*wn + 1j*wn*sqrt(1-zeta^2), ...
    -zeta*wn - 1j*wn*sqrt(1-zeta^2)];   % dominant pair
s_fast = [-12, -14];      % non-dominant: fast but not so fast that the
% gain (and control force) becomes excessive.
% Must be DISTINCT -- place() rejects repeats
% for a single-input system.

p_des = [s_dom, s_fast];
fprintf('\nDesired closed-loop poles:\n'); disp(p_des.');

% Desired characteristic polynomial (for the coefficient-matching slide)
chi_d = poly(p_des);
fprintf('Desired characteristic polynomial coefficients:\n'); disp(chi_d);

%% Compute the feedback gain
K_pp = place(A, B, p_des);
fprintf('K_pp = [%.2f  %.2f  %.2f  %.2f]\n', K_pp);

% Cross-check with Ackermann (valid here: single input)
K_ack = acker(A, B, p_des);
fprintf('K_ack = [%.2f  %.2f  %.2f  %.2f]   (max diff = %.2e)\n', ...
    K_ack, max(abs(K_pp - K_ack)));

% Verify the poles actually landed where we asked
fprintf('Achieved closed-loop poles:\n'); disp(eig(A - B*K_pp));

%% Closed-loop regulator simulation (same initial tilt as the open-loop test)
sys_cl_pp = ss(A - B*K_pp, B, C, D, ...
    'statename', states, 'inputname', inputs, 'outputname', outputs);

t  = 0:0.001:5;
x0 = [0; 0; 0.1; 0];                  % 0.1 rad initial tilt

[y_pp, t, x_pp] = initial(sys_cl_pp, x0, t);
u_pp = -K_pp * x_pp.';                % control effort

%% Performance metrics
tol   = 0.002;                        % 2% of the 0.1 rad initial condition
i_x   = find(abs(y_pp(:,1)) > tol, 1, 'last');
i_phi = find(abs(y_pp(:,2)) > tol, 1, 'last');
ts_pp = max(t(i_x), t(i_phi));

fprintf('\n--- Pole placement performance ---\n');
fprintf('Settling time      : %.2f s   (spec < 5 s)\n', ts_pp);
fprintf('Peak |phi|         : %.4f rad = %.2f deg   (spec < 0.35 rad)\n', ...
    max(abs(y_pp(:,2))), rad2deg(max(abs(y_pp(:,2)))));
fprintf('Peak |x|           : %.4f m   (spec < 0.5 m)\n', max(abs(y_pp(:,1))));
fprintf('Peak |u|           : %.2f N   (spec < 20 N)\n', max(abs(u_pp)));

%% Plots
figure;
subplot(3,1,1); grid on;
plot(t, y_pp(:,1), 'LineWidth', 1.4);
ylabel('x (m)'); title('Pole Placement: regulation from 0.1 rad tilt');
xlim([0 3]);

subplot(3,1,2); hold on; grid on;
plot(t, y_pp(:,2), 'r', 'LineWidth', 1.4);
yline(0.35, 'k:', 'LineWidth', 1.2);
yline(-0.35,'k:', 'LineWidth', 1.2);
ylabel('\phi (rad)'); xlim([0 3]);
legend('\phi', '\pm20\circ limit', 'Location', 'northeast');

subplot(3,1,3); grid on;
plot(t, u_pp, 'm', 'LineWidth', 1.4);
ylabel('u (N)'); xlabel('Time (s)'); xlim([0 3]);

%% Response characteristics
step_info = lsiminfo(y,t);
cart_info = step_info(1)   % characteristics for x    (use Max / MaxTime)
pend_info = step_info(2)   % characteristics for phi  (use Max / MaxTime)