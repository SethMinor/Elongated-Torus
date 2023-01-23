%% Numerical NLS vortex dynamics
% On the surface of the elongated torus
clear, clc;

% Fontsize, for plotting
fs = 14;

% NLS density
mu = 3;

% Grid size (NxN)
N = 80;

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
w1 = (0) + 1i*(4); % positive vortex
w2 = (0) + 1i*(-4); % negative vortex

% Complex flow potential
F =@(w,w1,w2) log(jacobitheta1((w-w1)./(2*c),p,cap)./jacobitheta1((w-w2)./(2*c),p,cap))...
    - (real(w1-w2)./(2*c^2*gr))*w;

% Phase of the flow field, Phi=Im(F)
phase =@(w,w1,w2) imag(F(w,w1,w2));

% Create a contour plot of the phase
figure (2)
Ugrid = linspace(-c*pi,c*pi,N);
Vgrid = linspace(c*gl,c*gr,N);
[Utemp,Vtemp] = meshgrid(Ugrid,Vgrid);
Z = phase(Utemp+1i*Vtemp,w1,w2);
contour(Utemp,Vtemp,Z,50)
colormap hsv;
axis equal;
title('Phase Contours')

% Create initial wave function
IC =@(w,w1,w2) sqrt(mu)*exp(1i*phase(w,w1,w2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w1)).^2 + (imag(w)-imag(w1)).^2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w2)).^2 + (imag(w)-imag(w2)).^2));

% Plot density of initial condition
figure (3)
psi_0 = IC(Utemp+1i*Vtemp,w1,w2);
Z = conj(psi_0).*psi_0;
surf(Utemp,Vtemp,Z)
shading interp;
colormap gray;
axis equal;
camlight
title('Initial Density')

% For plotting on el toroos:
% https://www.mathworks.com/help/matlab/visualize/representing-a-matrix-as-a-surface.html

%% Numerical integration
% RK-4 for time
dt = 0.07;
tf = 30;
N_time = floor(tf/dt);

% RHS parameters
D2 = Delta2d(N);

% Local scale factor
figure (4)
[Phi_temp,Theta_temp] = meshgrid(phi(Ugrid),theta(Vgrid)');
surf(Phi_temp,Theta_temp,lambda(Phi_temp,Theta_temp))
title('Local scale factor')

% Reshape initial condition array
seed = reshape(psi_0,[N^2,1]);

% Wrap lambda into vector
lambda_vec = reshape(lambda(Phi_temp,Theta_temp),[N^2,1]);

% Confirm IC
disp('Continue? Press ye olde key...')
pause;

% RK-4 for-loop
t = 0; % Initialize time
psi = seed; % Initialize wave function

for i = 0:N_time
    % Plot time step
    psi_temp = reshape(psi,[N,N]);
    density = conj(psi_temp).*psi_temp;
    figure (5)
    surf(Utemp,Vtemp,density)
    shading interp;
    colormap bone;
    axis equal;
    view(0,90)
    %camlight
    title("$t=$ "+t,'Interpreter','latex','FontSize',fs)
    xlim([-pi*c,pi*c])
    ylim([c*gl,c*gr])
    %pause(0.01)

    % Update using RK-4
    k1 = RHS(psi,D2,lambda_vec);
    k2 = RHS(psi + (dt/2)*k1,D2,lambda_vec);
    k3 = RHS(psi + (dt/2)*k2,D2,lambda_vec);
    k4 = RHS(psi + dt*k3,D2,lambda_vec);
    psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Update the time
    t = t + dt;
end


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
function F_of_psi = RHS(psi,D2,lambda_vec)
    F_of_psi = 1i*((D2*psi)./lambda_vec - (psi.^2).*conj(psi));
end
