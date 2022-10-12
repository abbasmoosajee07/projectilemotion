% Runge Kutta, 4th order, Function for projectile motion
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
function [rx_rk,ry_rk,vx_rk,vy_rk,tt_rk] = RungeKutta(v,theta,r_xy,dt)
rx_rk=r_xy(1,:); ry_rk=r_xy(2,:); % Assigning the specific X and Y coordinates from initial vector      
g=9.81;                  % Acceleration due to Gravity in m s^-2
vx_rk=v.*cosd(theta);    % X component of Velocity
vy_rk=v.*sind(theta);    % Y component of Velocity
t_max=(vy_rk+sqrt((vy_rk^2)+(2*g*ry_rk)))/g; % Maximum flight time
N=abs(t_max/dt); % Number of iterations loop should perform

    for ts = 1:N        %Runge Kutta code
        a1 = dt*vx_rk(ts);            b1 = dt*0;
        c1 = dt*(vy_rk(ts));          d1 = dt*(-g);
        a2 = dt*(vx_rk(ts) + b1/2);   b2 = dt*0;
        c2 = dt*(vy_rk(ts) + d1/2);   d2 = dt*(-g);
        a3 = dt*(vx_rk(ts) + b2/2);   b3 = dt*0;
        c3 = dt*(vy_rk(ts) + d2/2);   d3 = dt*(-g);
        a4 = dt*(vx_rk(ts) + b3);     b4 = dt*0;
        c4 = dt*(vy_rk(ts) + d3);     d4 = dt*(-g);
        
        rx_rk(ts+1) = rx_rk(ts)   + 1/6*(a1 + 2*a2 + 2*a3 + a4); % New X velocity calculated using respective integrated constants 
        vx_rk(ts+1) = vx_rk(ts)   + 1/6*(b1 + 2*b2 + 2*b3 + b4); % New X position calculated using respective integrated constants
        vy_rk(ts+1) = vy_rk(ts)   + 1/6*(d1 + 2*d2 + 2*d3 + d4); % New Y velocity calculated using respective integrated constants
        ry_rk(ts+1) = ry_rk(ts)   + 1/6*(c1 + 2*c2 + 2*c3 + c4); % New Y position calculated using respective integrated constants
        tt_rk(ts+1)=(ts+dt)*dt;   % Total array of time
    end
end
