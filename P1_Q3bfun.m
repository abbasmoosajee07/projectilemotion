% Labs and Data Analysis 2
% Portfolio 1
% Question 3: Part b
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
close all;clear;clc; %Clear Figures, Workspace, and Command Window
disp('Welcome to the Projectile Motion Plotter');
disp("Simulate a drag on ball of different materials in various fluids");
%% Constant Variables 
g=9.81;                 % Acceleration due to Gravity in m s^-2
theta=0;                % Angle reative to Horizontal has to be between 0-90 degrees
v=0;                    % Initial Velocity in m/s
r_xy=[0;100];            % Initial Position of ball as X and Y Coordinates
% Time Parameters
t_i=0;                  % Initial Time in seconds
dt=0.001;               % Time Step, the minimum time step which will accurately calculate terminal velocity of all particles is 0.035
t_max=1000;               % Max Time in seconds
T_list=(t_i:dt:t_max);
%% Defining characteristics for drag
% Particle Properties
P_mat=["Steel","Nylon"];% Particle Name
P_rho=[7750,1150];      % Density of particle material in kg/m^3
D=5E-03;                % Diameter of Sphere in m
P_a=pi*(D/2)^2;         % Projected Area of Sphere in m^2
Vol=4/3*(pi*(D/2)^3);   % Volume of Sphere in m^3
P_m=P_rho.*Vol;         % Mass of Particle in kg
Cd=0.47;                % Drag coefficient of a sphere ball is taken as 0.47
% Fluid Properties at RTP
F_nam=["Air";"Water";"Glycerine"]; % Names of the Fluids
F_rho=[1.25;1000;1250];            % Fluid Density in kg m^-3
%% Projectile Motion with Drag
tv_a  = sqrt((2*P_m.*g)./(F_rho.*P_a*Cd)); % Terminal Velocity of different balls calculated in various fluids theoretically
tv_wb=sqrt(((P_m*g)-(F_rho.*Vol*g))./(0.5*F_rho.*P_a*Cd)); % Terminal Velocity accounting for buyoancy

% [tt,r_y,r_x,v_y,a_y] = ProjectileDrag(P_m,P_a,F_rho,r_xy,dt,t_max)
[ts_sa,r_y_sa,r_x_sa,v_y_sa,a_y_sa] = ProjectileDrag(P_m(:,1),P_a,F_rho(1,:),r_xy,dt,t_max);  % Calling the Drag function for steel ball in air
[ts_sw,r_y_sw,r_x_sw,v_y_sw,a_y_sw] = ProjectileDrag(P_m(:,1),P_a,F_rho(2,:),r_xy,dt,t_max);  % Calling the Drag function for steel ball in water
[ts_sg,r_y_sg,r_x_sg,v_y_sg,a_y_sg] = ProjectileDrag(P_m(:,1),P_a,F_rho(3,:),r_xy,dt,t_max);  % Calling the Drag function for steel ball in glycerin

[ts_na,r_y_na,r_x_na,v_y_na,a_y_na] = ProjectileDrag(P_m(:,2),P_a,F_rho(1,:),r_xy,dt,t_max);  % Calling the Drag function for nylon ball in air
[ts_nw,r_y_nw,r_x_nw,v_y_nw,a_y_nw] = ProjectileDrag(P_m(:,2),P_a,F_rho(2,:),r_xy,dt,t_max);  % Calling the Drag function for nylon ball in water
[ts_ng,r_y_ng,r_x_ng,v_y_ng,a_y_ng] = ProjectileDrag(P_m(:,2),P_a,F_rho(3,:),r_xy,1E-3,5E-2); % Calling the Drag function for nylon ball in glycerin; dt and tmax are defined specifically
%% Displaying Projectile Flight Information
tv_sa = sprintf('The terminal velocity reached by steel ball in air %f m s^-1.\n',v_y_sa(:,end));       
tv_sw = sprintf('The terminal velocity reached by steel ball in water %f m s^-1.\n',v_y_sw(:,end));    
tv_sg = sprintf('The terminal velocity reached by steel ball in glycerin %f m s^-1.\n',v_y_sg(:,end));  
tv_na = sprintf('The terminal velocity reached by nylon ball in air %f m s^-1.\n',v_y_na(:,end));
tv_nw = sprintf('The terminal velocity reached by nylon ball in water %f m s^-1.\n',v_y_nw(:,end));
tv_ng = sprintf('The terminal velocity reached by nylon ball in glycerin %f m s^-1.\n',v_y_ng(:,end)); 
message=sprintf('%s',tv_sa,tv_sw,tv_sg,tv_na,tv_nw,tv_ng);disp(message);disp(tv_a);disp(tv_wb)

hFigure = figure;
% Figure Properties
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.0, 0.0, 01, 01]);
set(gcf, 'Name', 'Drag on objects in different fluids', 'NumberTitle', 'Off') 
subplot(3, 1, 1);
    hold on
    plot(ts_sa,v_y_sa,'.','Color',[0.4940 0.1840 0.5560]);
    plot(ts_na,v_y_na,'.','Color',[0.4660 0.6740 0.1880]);
    xlabel('Time (s)');  ylabel('Velocity (m s^-^1)')
    legend('Steel Ball','Nylon Ball','location','northeast')
    title('Drag for Balls in Air');
    xlim([0 10]);ylim([-35,0]);
    hold off
subplot(3, 1, 2);
    hold on
    plot(ts_sw,v_y_sw,'.','Color',[0.4940 0.1840 0.5560]);
    plot(ts_nw,v_y_nw,'.','Color',[0.4660 0.6740 0.1880]);
    xlabel('Time (s)');  ylabel('Velocity (m s^-^1)')
    title('Drag for Balls in Water');
    xlim([0,1.5]);ylim([-1.2,0.2]);
    hold off
subplot(3, 1, 3);
    hold on
    plot(ts_sg,v_y_sg,'.','Color',[0.4940 0.1840 0.5560]);
    plot(ts_ng,v_y_ng,'.','Color',[0.4660 0.6740 0.1880]);
    xlabel('Time (s)');  ylabel('Velocity (m s^-^1)');
    title('Drag for Balls in Glycerin');
    xlim([0,1.5]);ylim([-1.2,0.2]);
    hold off
