%% Finding fixed points on the elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;

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

%% Initial seeds for root-finding
% Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
% [-pi*c,pi*c] = [-33.2475, 33.2475]
% [cgl,cgr] = [-10.9479, 10.9479]
w1_0 = (3) + 1i*(0); % positive vortex
w2_0 = (-3) + 1i*(0); % negative vortex

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

% Real and imaginary parts of isothermal coords
u1_0 = real(w1_0);
v1_0 = imag(w1_0);

u2_0 = real(w2_0);
v2_0 = imag(w2_0);

% Initial seed
W = [u1_0, u2_0, v1_0, v2_0];

%% Root-finding
% Define the RHS
F =@(W) vortex_velocity_v2(0,[W(1), W(2), W(3), W(4)],0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
