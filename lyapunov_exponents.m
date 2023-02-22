%% Simulations on the elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;
%'Interpreter','latex','FontSize', fs

% Parameters
a = 12;
R = 11;
r = 3;
c = sqrt(R^2 - r^2);
myalpha = r/R;

%% Solving for isothermal coordinates
% Boundary condition
phi_0 = 0;

% Numerical integration using ode45
theta_span = [0, pi];
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);
[theta_raw,phi_raw] = ode15s(@(theta,phi) odefcn(theta,phi,a,R,r), theta_span, phi_0, options);

% Give these hosses odd symmetry 
theta_raw = [-flip(theta_raw(2:end)); theta_raw];
phi_raw = [-flip(phi_raw(2:end)); phi_raw];

% Quick plot to verify isothermal solutions look good
v0 =@(theta) 2*r*atan(sqrt((R-r)/(R+r))*tan(theta/2));

% Numerical isothermal coordinates
figure (1)
subplot(2,1,1)
plot(theta_raw, real(phi_raw))
hold on
plot(theta_raw, imag(phi_raw))
plot(theta_raw, -v0(theta_raw)/c,'--k')
xline(pi,'-')
xline(-pi,'-')
hold off
grid on

title('Numerical Isothermal Coordinates','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel('$\phi = f(\theta)$ solution','Interpreter','latex','FontSize', fs)
legend('Re$[f(\theta)]$','Im$[f(\theta)]$','$(-v_0)/c$','Interpreter','latex','FontSize', fs)

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp'); % phi = f(theta)
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Check derivative of numerical solution
D1 = differentiate(vhelper, theta_raw);
Df = fit(theta_raw, D1, 'cubicinterp');

subplot(2,1,2)
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
plot(theta_raw, abs(imag(-1i*r./gamma(phi_raw,theta_raw)) -  Df(theta_raw)),'.-')
grid on
set(gca, 'YScale', 'log')

title("Residual, max $=$ "+max(abs(imag(-1i*r./gamma(phi_raw,theta_raw)) -  Df(theta_raw))),'Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel("RHS - $f'(\theta)$",'Interpreter','latex','FontSize', fs)
legend('Im','Interpreter','latex','FontSize', fs)

% Conformal map
w =@(phi,theta) u(phi) + 1i*v(theta);

% Define gl and gr
gl = min(imag(phi_raw));
gr = max(imag(phi_raw));

% Defining the (u,v) to (phi,theta) map
ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
phi =@(u) u/c;
theta =@(v) ginv(-UVwrap(v, [-c*gr, c*gr])/c);

% Define the derivative of g_inverse
D2 = differentiate(ginv, imag(phi_raw));
Dginv = fit(imag(phi_raw), D2, 'cubicinterp');

%% Other parameters
% Gamma quantity
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);

% Local scale factor
lambda =@(phi,theta) gamma(phi,theta)./c;

% Jacobi theta function parameters
cap = 20;     % truncation error
p = exp(-gr); % nome (for periodicity)

%% Initial conditions for vortices
% Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
w1_0 = (-30) + 1i*(-10); % positive vortex
w2_0 = (6) + 1i*(5); % negative vortex

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

% real and imaginary parts of isothermal coords
u1_0 = real(w1_0);
v1_0 = imag(w1_0);

u2_0 = real(w2_0);
v2_0 = imag(w2_0);

%% Integrate the equations of motion
% Set total time and tolerances
t0 = 0;
tf = 500;
timespan = [t0, tf];
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

% Define the RHS
F =@(W) vortex_velocity_v2(0,[W(1), W(2), W(3), W(4)],0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);

% Numerical integration using ode45
y0 = [u1_0, u2_0, v1_0, v2_0]'; % ODE ICs
Q0 = eye(2*N); % Lyapunov, Q ICs
p0 = ones(2*N,1); % Lyapunov, rho_i ICs
Y0 = [y0; reshape(Q0,[(2*N)^2,1]); p0];

[t,y] = ode15s('vortex_velocity_lyapunov',timespan, Y0, options,...
    N, q, r, a, R, c, p, cap, theta, Dginv, gr, F);

%% Change coordinates from numerical solution
% Vortices in isothermal coordinates
U = y(:,1:N); % u-coords
V = y(:,(1+N):2*N); % v-coords

% Compute energy (before UVwrap)
[energy,classic,curve,quantum] = hamiltonian_v2(U,V,N,q,p,c,r,a,R,cap,phi,theta,gr);

% Unwrap
U = UVwrap(U, [-pi*c, pi*c]);
V = UVwrap(V, [c*gl, c*gr]);

% Full complex coordinate
W = U + 1i*V;

% Toroidal-poloidal coordinates
Phi = U./c;
Theta = [theta(V(:,1)), theta(V(:,2))];

% 3D Cartesian coordinates
X = (a+r*cos(Theta)).*cos(Phi);
Y = (R+r*cos(Theta)).*sin(Phi);
Z = r*sin(Theta);

%% Plot integrated solution
% Some nice RGB colores
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];

% 3D Cartesian surface plot of orbits
utorus = linspace(0,2*pi);
vtorus = linspace(0,2*pi);
[Utorus,Vtorus] = meshgrid(utorus,vtorus);
Xtorus = (a+r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R+r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);

figure (3)
surf(Xtorus, Ytorus, Ztorus,'FaceAlpha',0.3);
colormap('gray')
shading interp;
hold on
plot3(X,Y,Z)
hold off
grid on
axis equal
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)
title('3D Cartesian (Physical Path)','Interpreter','latex','FontSize',fs)

%% Display Hamiltonian
% Compute total energy
%[energy,classic,curve,quantum] = hamiltonian_v2(U,V,N,q,p,c,r,a,R,cap,phi,theta,gr);
energy_time = linspace(t0,tf,length(energy));

figure (4);

% Plot relative difference (E(t)-E(0))/E(0)
subplot(2,1,1)
plot(energy_time,((energy-energy(1))./energy(1)))
title('Energy of solution','Interpreter','latex')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\frac{H(t)-H_0}{H_0}$','Interpreter','latex')

% Plot each energy contribution
subplot(2,1,2)
plot(energy_time,energy,'-',energy_time,classic,'--')
hold on
plot(energy_time,curve,'--',energy_time,quantum,'--')
hold off
grid on
xlabel('$t$','Interpreter','latex')
ylabel('Energy contributions','Interpreter','latex')
legend('Total, $H$','Classic','Curvature','Quantum','Interpreter','latex')

%% Compute Lyapunov exponents
% Roll through 'y_n' solution and store {rho_i} values
% [rho_1, ..., rho_4] list = y(:, 1 + 2N + (2N)^2 : end)
rho_list = y(:, 1+2*N+(2*N)^2 : end);

% Compute the full spectrum
N_lyapunov = 200; % number to back-average
Lyapunov = zeros(1,2*N);
for i = 1:2*N
    Lyapunov(i) = mean(rho_list(end-N_lyapunov:end,i)./t(end-N_lyapunov:end));
end

%% Plot the Lyapunov spectrum
plot_every = 5; % Speed-up for plotting

figure (5)
sgtitle('Lyapunov Spectrum','Interpreter','latex','FontSize',fs+2)

subplot(4,1,1)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,1)./t(1:plot_every:end))
hold on
yline(Lyapunov(1),'--r')
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{1}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_1=$ "+Lyapunov(1),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Interpreter','latex')
ylim([Lyapunov(1)-0.5, Lyapunov(1)+0.5])

subplot(4,1,2)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,2)./t(1:plot_every:end))
hold on
yline(Lyapunov(2),'--r')
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{2}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_2=$ "+Lyapunov(2),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Interpreter','latex')
ylim([Lyapunov(2)-0.5, Lyapunov(2)+0.5])

subplot(4,1,3)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,3)./t(1:plot_every:end))
hold on
yline(Lyapunov(3),'--r')
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{3}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_3=$ "+Lyapunov(3),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Interpreter','latex')
ylim([Lyapunov(3)-0.5, Lyapunov(3)+0.5])

subplot(4,1,4)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,4)./t(1:plot_every:end))
hold on
yline(Lyapunov(4),'--r')
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{4}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_4=$ "+Lyapunov(4),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Interpreter','latex')
ylim([Lyapunov(4)-0.5, Lyapunov(4)+0.5])

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
