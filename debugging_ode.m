%% Simulations on the elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;
%'Interpreter','latex','FontSize', fs

% Parameters
% SEEMS LIKE THIS CODE MATCHES OG CODE BETTER FOR LARGER ALPHA
% BATTLE BETWEEN ISOTHERMAL AND THIS (?)
% Symplectic integrator as a possible fix (?)
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
subplot(3,1,1)
plot(theta_raw, real(phi_raw))
hold on
plot(theta_raw, imag(phi_raw))
plot(theta_raw, -v0(theta_raw)/c,'--k')
hold off
grid on

title('Isothermal Coordinates ODE','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel('$\phi = f(\theta)$ solution','Interpreter','latex','FontSize', fs)
legend('Re$[f(\theta)]$','Im$[f(\theta)]$','$(-v_0)/c$','Interpreter','latex','FontSize', fs)

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp');
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Check residuals/derivatives
D1 = differentiate(vhelper, theta_raw);
Df = fit(theta_raw, D1, 'cubicinterp');

subplot(3,1,2)
% RHS = -i*r/gamma
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
plot(theta_raw, real(-1i*r./gamma(phi_raw,theta_raw)),'--')
hold on
plot(theta_raw, imag(-1i*r./gamma(phi_raw,theta_raw)),'--')
plot(theta_raw(2:end), real(diff(phi_raw)./(2*pi/length(theta_raw))))
plot(theta_raw, Df(theta_raw))
grid on
% Circular torus comparison
dv0 =@(theta) -r./(R + r*cos(theta));
plot(theta_raw, dv0(theta_raw), '--k')
hold off
% LOOKS LIKE Im(f) = -r/gamma (?!)

title('Derivative Matching','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel("Im$[f'(\theta)]$, solution derivative",'Interpreter','latex','FontSize', fs)
legend('RHS (Re)','RHS (Im)','Numerical (Re)', 'Numerical (Im)',"$(-v_0')/c$",'Interpreter','latex','FontSize', fs)

subplot(3,1,3)
title('Residual','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel("Im$[f'(\theta)]$, solution derivative",'Interpreter','latex','FontSize', fs)
legend('RHS (Re)','RHS (Im)','Numerical (Re)', 'Numerical (Im)',"$(-v_0')/c$",'Interpreter','latex','FontSize', fs)

% Conformal map
w =@(phi,theta) u(phi) + 1i*v(theta);

% Defining the (u,v) to (phi,theta) map
ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
phi =@(u) u/c;
theta =@(v) ginv(-v/c);

% Define the derivative of g_inverse
D2 = differentiate(ginv, imag(phi_raw));
Dginv = fit(imag(phi_raw), D2, 'cubicinterp');

% Define gl and gr
gl = min(imag(phi_raw));
gr = max(imag(phi_raw));

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
