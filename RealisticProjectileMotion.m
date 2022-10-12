% Labs and Data Analysis 2
% Portfolio 1: Combined(most realistic projectile motion)
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
close all;clear;clc; %Clear Figures, Workspace, and Command Window
disp('Welcome to the Projectile Motion Plotter');
%% User Defined Data
Planet= ['Earth=1', 'Mars=2','Jupiter=3','Sun=4'];     % Select Planet on whcih projectile is launched
Gravity=[9.81; 3.27; 24.87; 274];  g=Gravity(1,:);     % Input corresponding number from Planet array into g where all values are in m s^-2

P_nam = ['Steel=1', 'Nylon=2','Leather=3'];            % Select particle of whcih density is made off
P_Rho = [7500; 1250; 802];P_rho=P_Rho(2,:);            % Density of particle material in kg/m^3

F_nam=["Air";"Water";"Glycerine";"Random Fluid"];      % Names of the Fluids
F_Rhom=[1.25;1000;1250;2000];  F_rho=F_Rhom(1,:);      % Fluid Density in kg m^-3

r_x=(0);                            % Select initial X position
r_y=(100);                          % Select initial Y position

% User Input Angle
    theta=(30);   % Angle reative to Horizontal can be changed, given its between 0-90degrees
    if (0 <= theta && theta <=90 ) % If condition that evaluates if angle is between 0 and 90 degrees
    else
        fprintf('\n Angle input was not within range. \n A default of 45 degrees will be used instead \n ');
        pause(2),theta=45; 
    end

% User Input Velocity
v=(40);       % Initial Velocity in m/s max y velocity of 28.875
    if (0 <= v && v <=(28.875/sind(theta)))
    else
        fprintf('\n Velocity input was not within valid ranger. \n A default of 10 m/s will be used instead \n ');
        pause(2),v=10; 
    end
    
D=(0.15);                % Diameter of Sphere in m
dt=0.001;               % Step size
t_max=1000000;          % Maximum time span of calculations

%% Type of particle  motion
fprintf("Select type of particle motion. \n   For Simple projectile enter:    1 " + ...
    " \n   For a bouncing particle enter:  2 \n   For Drag enter: 3 \n   For Drag with bounce enter: 4 \n")
choice_pm=input('Input Projectile motion type:  ');
    if choice_pm == 1
        Cd=0.00;                % Both Drag and impact collision are negated and so 
        Re=0.00;                % only projectile motion is considered
    elseif choice_pm == 2
        Cd=0.00;                % This condition will calculate the motion for a bouncing 
        Re=0.59;                % particle
    elseif choice_pm == 3
        Cd=0.47;                % This condition will calculate the projectile motion for
        Re=-1.00;               % a particle with drag and so only Cd is take
    elseif choice_pm == 4
        Cd=0.47;                % Drag coefficient of a sphere ball is taken as 0.47
        Re=0.59;                % Restitution co-efficient, non-dimensional
    elseif choice_pm > 4
        fprintf("Please input a valid number"), pause(2)
        run('RealisticProjectileMotion.m')
    end

%% Initial Calculations
P_a=pi*(D/2)^2;         % Projected Area of Sphere in m^2
Vol=4/3*(pi*(D/2)^3);   % Volume of Sphere in m^3
P_m=P_rho.*Vol;         % Mass of Particle in kg
v_x=v*cosd(theta);      % X component of Velocity (m/s)
v_y=v*sind(theta);      % Y component of Velocity (m/s)
Fb=F_rho.*g.*Vol;       % Buoyancy force  (N)
Fw=P_m*g;               % Particle Weight (N)
FD_x = 0.5*Cd*F_rho*P_a*v_x^2;    % Inital Drag force in X
FD_y = 0.5*Cd*F_rho*P_a*v_y^2;    % Inital Drag force in Y
TFr_x=FD_x; TFr_y=Fw-(FD_y+Fb);   % Resultant Force in X and Y
a_x = TFr_x/P_m; a_y = TFr_y/P_m; % Inital Acceleration in X and Y
GPE= P_m*g*r_y;KE = P_m*0.5*v_y.^2;TE = GPE+KE;   % Initial Energy of Partile
N=abs(t_max/dt);        % Calculating the Number of iterations loop needs to perform   

%% Main Loop for projectile motion
for ts=1:N
    ts=ts+1;

    v_x(ts) = v_x(ts-1)+(-a_x(ts-1))*dt;  % New X velocity, updated using x acceleration
    v_y(ts) = v_y(ts-1)+(-a_y(ts-1))*dt;  % New Y velocity, updated using y acceleration
    
    r_x(ts) = r_x(ts-1)+v_x(ts-1)*dt;     % New X position, updated using x velocity
    r_y(ts) = r_y(ts-1)+v_y(ts-1)*dt;     % New Y position, updated using y velocity

    GPE(ts) = P_m*g*r_y(ts-1);            % Gravitational Potential Energy calculated using m*g*h in J
    KE(ts)  = P_m*0.5*v_y(ts-1)^2;        % Kinetic Energy calculated using 0.5*m*v^2 in J
    TE(ts)  = GPE(ts-1)+KE(ts-1);         % Total Energy of Particle, sum of GPE and KE
    
    FD_x(ts) = 0.5*Cd*F_rho*P_a*v_x(ts)^2;% 
    FD_y(ts) = 0.5*Cd*F_rho*P_a*v_y(ts)^2;% 

    TFr_x(ts)=FD_x(ts); a_x(ts) = TFr_x(ts)/P_m;
    TFr_y(ts)=Fw-(FD_y(ts)+Fb); a_y(ts) = TFr_y(ts)/P_m;
    
        if r_y(ts-1) <0                % If condition that determines when projectile reaches the ground
            v_y(ts) =-Re* v_y(ts-1);   % Find new y velocity, as particle bounces back
            r_y(ts) = 0;               % The  y position, as particle bounces back is 0
        else
            v_y(ts) = v_y(ts-1)+(-a_y(ts-1))*dt;  % Find new y velocity, if no impact occurs
        end                            
     tt(ts)=((ts-1)+dt)*dt;
          if TE(ts) < P_m*g            % If condition determines if total energy of particle is less than 0.05 Energy required to leave sureface collison is factor of mass, weight and gravity
             break                     % IF true, for loop is broken
          else
             if v_y(ts)==v_y(ts-1) % If condition that determnes when two consecutive velocities in the
                 break             % are the same, and so terminal thus breaking the for loop
             end
         end                           
end

%% Displaying Projectile Flight Information
hFigure = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.0, 0.0, 01, 01]);
set(gcf, 'Name', 'Projectile motion for a bouncing ball with drag', 'NumberTitle', 'Off') 
subplot(3, 3, 1);
    plot(r_x,r_y,'.')
    hold on 
    title('Height against horizontal displacement')
    xlabel('Horizontal Distance(m)'); ylabel('Height(m)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(3, 3, 2);
    plot(tt,r_x,'.')
    hold on 
    title('Horizontal Displacement against time')
    xlabel('Time(s)'); ylabel('Horizontal Displacement(m)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(3, 3, 3);
    plot(tt,r_y,'.')
    hold on 
    title('Vertical height against time')
    xlabel('Time(s)'); ylabel('Height(m)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(3, 2, 3);
    plot(tt,v_y,'.')
    hold on 
    title("Velocity against Time"); 
    xlabel('Time(s)'); ylabel('Velocity(m s^-^1)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(3, 2, 4);
    hold on
    plot(tt,a_y,'k.-');
    title('Acceleration against time')
    xlabel('Time (s)');ylabel('Acceleration (m s^-^2)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off
subplot(3, 1, 3);
    hold on 
    plot(tt,GPE,'.','Color',[0.4660, 0.6740, 0.1880]);
    plot(tt,KE, '.','Color',[0.8500, 0.3250, 0.0980]);
    plot(tt,TE,'.-','Color',[0 0 0])
    title("Energy of Particle through out its motion");
    legend("GPE",'KE','TE','location','northeast')
    xlabel('Time(s)'); ylabel('Energy(J)');
    g = gca; g.XAxisLocation = 'origin'; g.YAxisLocation = 'origin';
    hold off

m1 = sprintf('\n      The total time of travel is %f seconds.\n', tt(:,end));
m2 = sprintf('      The total distance travelled by particle is %f meters.\n',r_x(:,end));
m3 = sprintf('      The maximum height reached by particle is %f meters.\n',max(r_y)); 
message = sprintf('%s', m1, m2,m3); disp(message)
