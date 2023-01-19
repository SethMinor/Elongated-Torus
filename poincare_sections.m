% Poincare sections for elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;

% Parameters
a = 13;
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

figure (1)
plot(theta_raw, real(phi_raw))
hold on
plot(theta_raw, imag(phi_raw))
plot(theta_raw, -v0(theta_raw)/c,'--k')
hold off
grid on

title('Isothermal Coordinates ODE','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$','Interpreter','latex','FontSize', fs)
ylabel('$\phi = f(\theta)$ solution','Interpreter','latex','FontSize', fs)
legend('Re$[f(\theta)]$','Im$[f(\theta)]$','$-v_0/c$','Interpreter','latex','FontSize', fs)

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp');
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Conformal map
w =@(phi,theta) u(phi) + 1i*v(theta);

% Defining the (u,v) to (phi,theta) map
ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
phi =@(u) u/c;
theta =@(v) ginv(-v/c);

% Define the derivative of g_inverse
D = differentiate(ginv, imag(phi_raw));
Dginv = fit(imag(phi_raw), D, 'cubicinterp');

% Define gl and gr
gl = min(imag(phi_raw));
gr = max(imag(phi_raw));

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
% [-pi*c,pi*c] = [-33.2475, 33.2475]
% [cgl,cgr] = [-10.5830, 10.5830]
w1_0 = (5.0) + 1i*(3.0); % positive vortex
w2_0 = (4.0) + 1i*(-2.0); % negative vortex

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

%% Compute Poincare sections
% ode45 with events function
t0 = 0;
tf = 1000;
timespan = [t0, tf];
%options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', @EventsFcn);

% Numerical integration using ode45 or ode15s
y0 = [u1_0, u2_0, v1_0, v2_0];
%[t,y] = ode15s('vortex_velocity_v2',timespan, y0, options, N, q, r, a, R, c, p, cap, theta, Dginv, gr);
F =@(y) vortex_velocity_v2(0,y,0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);
%[T,Y,Te,Ye,Ie]=ode45(@(t,y) LorenzEqs(t,y,params),[0,tpoincare],IC,options2);

[T,Y,Te,Ye,Ie] = ode15s(@(t,y) vortex_velocity_v2(0,y,0,N,q,r,a,R,c,p,cap,theta,Dginv,gr),timespan, y0, options);

%% Change coordinates from numerical solution
% Vortices in isothermal coordinates
U = Y(:,1:N); % u-coords
V = Y(:,(1+N):2*N); % v-coords
W = U + 1i*V;
% mod(...-pi, pi) (?)

% Toroidal-poloidal coordinates
Phi = U./c;
Phi = unwrap(Phi); % reduce mod 2pi

Theta = [theta(V(:,1)), theta(V(:,2))];
Theta = unwrap(Theta);

%% Plot integrated solution
% Some nice RGB colores
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];

figure (2)

% Isothermal orbit
subplot(2,1,1)
plot(U,V)
grid on
xlabel('$u = $Re$(w)$','Interpreter','latex','FontSize',fs)
ylabel('$v = $Im$(w)$','Interpreter','latex','FontSize',fs)
title('Isothermal Coordinates','Interpreter','latex','FontSize',fs)

% Toriodal-poloidal orbit
subplot(2,1,2)
plotwrapped(Phi(:,1),Theta(:,1),1, [-pi pi],[-pi pi], 0.05, bluey)
hold on
plotwrapped(Phi(:,2),Theta(:,2),1, [-pi pi],[-pi pi], 0.05, orangu)
hold off
xlabel('$\phi$','Interpreter','latex','FontSize',fs)
ylabel('$\theta$','Interpreter','latex','FontSize',fs)
title('Toroidal-Poloidal Coordinates','Interpreter','latex','FontSize',fs)
xlim([-pi pi])

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Poincare section slice
% UPDATE THIS BL WHEN R,r CHANGE
function [position,isterminal,direction] = EventsFcn(~,y)
  U = y(1:2);
  %V = y(3:4);
  R = 11;
  r = 3;
  c = sqrt(R^2-r^2);
  Phi = U./c;
  Phi = unwrap(Phi);
  % Theta = ...
  position = Phi(2) - pi;
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end
