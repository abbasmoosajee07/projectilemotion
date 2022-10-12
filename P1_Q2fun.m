% Labs and Data Analysis 2
% Portfolio 1
% Question 2
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
close all;clear;clc  %Clear Figures, Workspace, and Command Window
disp('Welcome to the Projectile Motion Plotter')
disp("Compare the motion of projectiles launched at 30 and 60 degrees.");
%% Dependent Variables in Projectile motion
theta=[60;30];      % Angle reative to Horizontal in degrees.
r_xy=[0;0];         % Initial Position of Particle %Ground Level coordinates
v=10;               % Initial velocity of 10 m s^-1
dt=0.01;         % Time Step
%% Calling the functions
% [r_x,r_y,v_x,v_y,r_xa,r_ya,tt] = Eulersfunction(v,theta,r_xy,dt)
% [rx_rk,ry_rk,vx_rk,vy_rk,tt_rk] = RungeKutta(v,theta,r_xy,dt)
[r_xe60,r_ye60,v_xe60,v_ye60,~,~,t_me60] = Eulersfunction(v,theta(1,:),r_xy,dt);  % Calling Euler's for 60 degrees
[rx_rk_60,ry_rk_60,vx_rk_60,vy_rk_60,tt_rk60] = RungeKutta(v,theta(1,:),r_xy,dt); % Calling Runge Kutta for 60 degrees

[r_xe30,r_ye30,v_xe30,v_ye30,~,~,t_me30] = Eulersfunction(v,theta(2,:),r_xy,dt);  % Calling Euler's for 30 degrees
[rx_rk_30,ry_rk_30,vx_rk_30,vy_rk_30,tt_rk30] = RungeKutta(v,theta(2,:),r_xy,dt); % Calling Runge Kutta for 30 degrees
%% Displaying Projectile Flight Information
hFigure = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.0, 0.0, 0.5, 0.5]);
set(gcf, 'Name', 'Projectile motion of 30 and 60 degrees', 'NumberTitle', 'Off') 
subplot(2, 1, 1);
plot(r_xe60,r_ye60,'r.')
hold on    
    plot(r_xe30,r_ye30,'m.')
    plot(rx_rk_60,ry_rk_60,'k-',rx_rk_30,ry_rk_30,'b-')
    title("Height against displacement of particles launched at different angles");
    legend("60 degree Euler's Method", "30 degree Euler's Method", ...
           "60 degree Runge Kutta's Method","30 degree Runge Kutta Method",'location','northeast')
    xlabel('Distance(m)'); ylabel('Height(m)'); grid on
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(2, 1, 2);
hold on  
plot(t_me60,r_ye60,'r.')
    plot(t_me30,r_ye30,'m.')
    plot(tt_rk60,ry_rk_60,'k-',tt_rk30,ry_rk_30,'b-')
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    xlabel('Time(s)'); ylabel('Height(m)'); grid on
hold off
