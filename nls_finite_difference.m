%% Numerical NLS vortex dynamics
% On the surface of the elongated torus
clear, clc;

% Fontsize, for plotting
fs = 14;

% NLS density
mu = 1;

% Torus parameters
a = 11;
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

%% Initial conditions for wave function
% Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
w1 = (12) + 1i*(2); % positive vortex
w2 = (-12) + 1i*(-5); % negative vortex

% Complex flow potential
F =@(w,w1,w2) log(jacobitheta1((w-w1)./(2*c),p,cap)./jacobitheta1((w-w2)./(2*c),p,cap))...
    - (real(w1-w2)./(2*c^2*gr))*w;

% Phase of the flow field, Phi=Im(F)
phase =@(w,w1,w2) imag(F(w,w1,w2));

% Create a contour plot of the phase
figure (2)
x = linspace(-c*pi,c*pi,500);
y = linspace(c*gl,c*gr,500);
[X,Y] = meshgrid(x,y);
Z = phase(X+1i*Y,w1,w2);
contour(X,Y,Z,50)

% Create initial wave function
psi_0 =@(w,w1,w2) sqrt(mu)*exp(1i*phase(w,w1,w2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w1)).^2 + (imag(w)-imag(w1)).^2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w2)).^2 + (imag(w)-imag(w2)).^2));

% Plot density of initial condition
figure (3)
Z = psi_0(X+1i*Y,w1,w2);
Z = conj(Z).*Z;
surf(X,Y,Z)
shading interp;
colormap copper;
axis equal;




%% Helper functions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Wrap U,V to interval
function wrapped = UVwrap(array, interval)
    wrapped = mod(array - interval(1), range(interval)) + interval(1);
end

% Creates Delta1D matrices
function D1 = Delta1d(N)
    D1 = diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    D1(end,1) = 1;
    D1(1,end) = 1;
end

% Places the Iu and Id matrices in Delta2d
function PlacedI = PlaceI(N)
    % Create the upper I-matrix portion
    IuBlock = {0};
    for m = 1:N
        IuBlock{m} = eye(N);
    end
    IuBlock = blkdiag(IuBlock{1:end});
    IuBlock = circshift(IuBlock,N,2);

    % Create the lower I-matrix portion
    IdBlock = {0};
    for m = 1:N
        IdBlock{m} = eye(N);
    end
    IdBlock = blkdiag(IdBlock{1:end});
    IdBlock = circshift(IdBlock,N,2)';

    % Return the combined blocks
    PlacedI = IuBlock + IdBlock;
end

% Returns the Delta2D matrix
function D2 = Delta2d(N)
    % Add the Delta1D blocks that run down the diagonal
    Block = {0};
    for m = 1:N
        Block{m} = Delta1d(N);
    end
    D2 = blkdiag(Block{1:end});

    % Add the Iu and Id blocks on the off-diagonals
    PlacedI = PlaceI(N);

    % Deliver the final spicy meatball
    D2 = D2 + PlacedI;
end

% Return the RHS
% LAMBDA(PSI) OR LAMBDA(UN,VN) KINDA THING (?)
function F_of_psi = RHS(psi,D2,lambda)
    F_of_psi = 1i*((D2*psi)./lambda(real(psi),imag(psi)) - (psi.^2).*conj(psi));
end

% Convert a numerical solution vector into a 2D grid
% function wavefunction = psi2D(v,Nphi,Ntheta)
%     wavefunction = reshape(v,[Nphi,Ntheta]);
% end
