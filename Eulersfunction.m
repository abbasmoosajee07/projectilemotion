% Euler's Function for projectile motion
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
function [r_x,r_y,v_x,v_y,r_xa,r_ya,tt] = Eulersfunction(v,theta,r_xy,dt)
r_x=r_xy(1,:); r_y=r_xy(2,:); % Assigning the specific X and Y coordinates from initial vector
g=9.81;                       % Acceleration due to Gravity in m s^-2
v_x=v.*cosd(theta);           % X component of Velocity
v_y=v.*sind(theta);           % Y component of Velocity
t_max=(v_y+sqrt((v_y^2)+(2*g*r_y)))/g; % Maximum time of flight accounting for intial particle height
N=abs(t_max/dt);              % Calculating the Number of iterations loop needs to perform
%% Analytic Solution
t_i=0;T_list=(t_i:dt:t_max);              % Complete array of times with an initial of zero
r_ya = r_y+(v_y.*T_list-(g.*T_list.^2/2));% Analytical Y Position (m)
r_xa = r_x+(v_x.*T_list);                 % Analytical X Position (m)
%% Euler's Numerical method
for ts = 1:N %
    ts=ts+1;
%   v_x(ts) = v_x(ts-1);               % New X Velocity, but no air resistance
    r_x(ts) = r_x(ts-1) + dt*v_x;      % New X Position (m)
    v_y(ts) = v_y(ts-1) + dt*(-g);     % New Y Velocity (m/s)
    r_y(ts) = r_y(ts-1) + dt*v_y(ts);  % New Y Position (m)
    tt(ts)  = (ts-1+dt)*dt;            % Total time of flight (s)
end                                    % End for loop
end                                    % End Function