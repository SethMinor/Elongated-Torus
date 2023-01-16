function dydt = vortex_velocity_v2(~,y0,~,N,q,r,a,R,c,p,cap,theta,Dginv,gr)
%% This function computes the physical velocity of a set of vortices,

% INPUTS
% y0 = initial conditions vector in isothermal, y0 = [u1_0 u2_0 v1_0 v2_0]
% N = number of vortices
% q = vector of vortex charges (should sum to 0)
% r,a,R,c = geometric parameters for torus
% p = nome of first Jacobi theta function (should satisfy |p| < 1)
% cap = truncation term of sum in first Jacobi theta function
% phi = phi(u) map
% theta = theta(v) map

% OUTPUTS
% omega = Vy + iVx = complex physical velocity of vortex cores
% dydt = right-hand side of ODE

%% Right-hand side of ODE
% Real and imaginary parts of isothermal vortex coords
y0 = y0(:);
U = y0(1:N);
V = y0((N+1):2*N);
W = U + 1i*V;

% Gamma and derivatives
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
gamma_phi = @(phi,theta) (a-R)*(a+R+2*r*cos(theta)).*cos(phi).*sin(phi)./gamma(phi,theta);
gamma_theta =@(phi,theta) -r*(a+R+2*r*cos(theta)+(R-a)*cos(2*phi)).*sin(theta)./(2.*gamma(phi,theta));

% Create Phi=Phi(U) and Theta=Theta(V)
Phi = U./c;
Theta = theta(V);

% Local scale factor terms
L = gamma(Phi,Theta)/c;
L_u = gamma_phi(Phi,Theta)/(c^2);
L_v = -(gamma_theta(Phi,Theta)/(c^2)).*Dginv(-V/c);

% Initialize omega vector
omega = zeros(N,1);

% Loop over each vortex
for n=1:N
    % Add curvature and quantum contributions
    curve = (q(n)./(2*L(n))).*(1i*L_v(n) - L_u(n));
    quantum = -(1/(2*c^2*gr))*(q(n).*U(n));
    omega(n) = omega(n) + curve + quantum;

    % Add induced velocity contributions (from other vortices)
    for m = 1:N
        if (m ~= n)
            T = jacobitheta1((W(n)-W(m))./(2*c),p,cap); % jacobi theta
            Tprime = Djacobitheta1((W(n)-W(m))./(2*c),p,cap);
            fmn = 1/(2*c)*(Tprime./T)-(1/(2*c^2*gr))*U(m);
            omega(n) = omega(n) + q(m).*fmn;
        end
    end

    omega(n) = (1./L(n)).*omega(n);
end

vdots = real(omega)./L;
udots = imag(omega)./L;

dydt = [udots; vdots];
end