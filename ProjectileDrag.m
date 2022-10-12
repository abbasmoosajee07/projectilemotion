% Function for Drag
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
function [tt,r_y,r_x,v_y,a_y] = ProjectileDrag(P_m,P_a,F_rho,r_xy,dt,t_max)
r_x=r_xy(1,:); r_y=r_xy(2,:); % Assigning the specific X and Y coordinates from initial vector
g=9.81;               % Acceleration due to Gravity in m s^-2
theta=0;              % Angle reative to Horizontal has to be between 0-90degrees
v=0;                  % Initial Velocity in m/s
v_x=v*cosd(theta);    % X component of Velocity (m/s)
v_y=v*sind(theta);    % Y component of Velocity (m/s)
Cd=0.47;              % Drag coefficient of a sphere is 0.47
R=sqrt(P_a/pi); Vol_dp=4/3*(pi*(R)^3); % Volume displaced m^3
Fb=F_rho.*g.*Vol_dp;  % Buoyancy force  (N)
Fw=P_m*g;             % Particle Weight (N)
N=abs(t_max/dt);      % Calculating the Number of iterations loop needs to perform   

%%Initial Force Calculations
FD_x = 0.5*Cd*F_rho*P_a*v_x^2;      % Inital Drag force in X
FD_y = 0.5*Cd*F_rho*P_a*v_y^2;      % Inital Drag force in Y
TFr_x=FD_x; TFr_y=Fw-(FD_y+Fb);     % Resultant Force in X and Y
a_x = TFr_x/P_m; a_y = TFr_y/P_m;   % Inital Acceleration in X and Y

for ts=1:N
    ts=ts+1;  % Begins the loop with initial ts value

    v_x(ts) = v_x(ts-1)+(-a_x(ts-1))*dt;  % New X velocity
    v_y(ts) = v_y(ts-1)+(-a_y(ts-1))*dt;  % New Y velocity
    
    r_x(ts) = r_x(ts-1)+v_x(ts-1)*dt;     % New X position
    r_y(ts) = r_y(ts-1)+v_y(ts-1)*dt;     % New Y position  

    FD_x =  0.5*Cd*F_rho*P_a*v_x(ts)^2;   % Drag force in X
    FD_y =  0.5*Cd*F_rho*P_a*v_y(ts)^2;   % Drag Force in Y

    TFr_x=FD_x;  a_x(ts) = TFr_x/P_m;        % Total force and acceleration in X
    TFr_y=Fw-(FD_y+Fb); a_y(ts) = TFr_y/P_m; % Total force and acceleration in Y, used arrays so that its valid for a projectile

    tt(ts)=((ts-1)+dt)*dt;                   % Array of Time Values
end
end
