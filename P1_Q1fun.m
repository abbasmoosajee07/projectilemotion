% Labs and Data Analysis 2
% Portfolio 1
% Question 1
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
close all;clear;clc %Clear Figures, Workspace, and Command Window
disp('Welcome to the Projectile Motion Plotter')
disp("A comparison of Projectile Motion using Euler's and Analytical Methods."); 
%% Define Initial Particle Data
r_x=(0);% Define initial position as X-Coordinate        
r_y=(0); % Define initial position as  Y-Coordinate
r_xy=abs([r_x;r_y]); % Converts coordinates to positive

% User Input Angle
theta=(60);   % Angle reative to Horizontal can be changed, given its between 0-90degrees
if (0 <= theta && theta <=90 ) % If condition that evaluates if angle is between 0 and 90 degrees
else
    fprintf('\n Angle input was not within range. \n A default of 45 degrees will be used instead \n ');
    pause(2), theta=45; 
end

% User Input Velocity
v=(10);       % Initial Velocity in m/s max
if (0 <= v )
else
    fprintf('\n Velocity input was not within valid ranger. \n A default of 10 m/s will be used instead \n ');
    pause(2),v=10; 
end
dt=0.01;               % Step size

%% Time Parameters
g=9.81;    % Earths Gravity in m s^-2
t_i=0;     % Initial Time
[r_x,r_y,~,~,r_xa,r_ya,tt] = Eulersfunction(v,theta,r_xy,dt); % Calling Eulers function
%% Displaying Projectile Flight Information
hfigure=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.0, 0.0, 01, 01]);
set(gcf, 'Name', "Projectile motion using Euler's", 'NumberTitle', 'Off') 
subplot(1,1,1)
hold on    
    plot(r_x,r_y,'ro')
    plot(r_xa,r_ya,'k-')
    title("Projectile Motion: Comparing Euler's with Analytical");
    legend("Euler's Method","Analytical Method",'location','northwest')
    xlabel('Distance(m)'); ylabel('Height(m)');
    ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
hold off

m1 = sprintf('\n      The total time of travel is %f seconds.\n', tt(:,end));
m2 = sprintf('      The total distance travelled by particle is %f meters.\n',r_x(:,end));
m3 = sprintf('      The maximum height reached by particle is %f meters.\n',max(r_y)); 
message = sprintf('%s', m1, m2,m3); disp(message)
