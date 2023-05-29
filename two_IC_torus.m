%% ODE for integrating two ICs simultaneously

function [dYdt] = two_IC_torus(~,y,params)
% Name parameters
[N,q,p,c,r,a,R,cap,~,theta,gr,Dginv] = params{:};

y1 = y(1:4); % First IC
y2 = y(5:8); % Second IC

dYdt(1:4) = vortex_velocity_v2(0,y1,0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);
dYdt(5:8) = vortex_velocity_v2(0,y2,0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);
dYdt = dYdt';
end
