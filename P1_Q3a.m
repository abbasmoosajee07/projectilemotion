% Labs and Data Analysis 2
% Portfolio 1
% Question 3: Part a
%
% Written by:   Abbas Moosajee
% Written on:   11/03/2022
% Contact:      AHM080@student.bham.ac.uk
close all;clear; %Clear Figures, Workspace, and Command Window
disp('Welcome to the Projectile Motion Plotter')
disp("Where the Euler's and Analytical Methods are used to compare Projectile motion."); 
%% User Defined Data
r_x=(0);% Define initial position as X-Coordinate        
r_y=(250); % Define initial position as  Y-Coordinate
dt=0.5;                % Time step
t_max=1E+12;            % Maximum time over which moion is simulated
N=abs(t_max/dt);        % Maximum number of iterations that for loop must perform
g=9.81;                 % Acceleration due to Gravity in m s^-2 when particle is dropped on earth

%% Projectile Variables
Re=0.59;                % Restitution co-efficient, non-dimensional
theta=0;                % Angle reative to Horizontal is 0 degrees
v=0;                    % As particle is dropped, inital velocity is zero
v_x=v.*cosd(theta);     % X component of Velocity in m s^-1
v_y=v.*sind(theta);     % Y component of Velocity in m s^-1
m=0.15;                 % Mass of hard steel ball in kg
GPE= m*g*r_y;KE = m*0.5*v_y.^2;TE = GPE+KE;   % Initial Energy of Partile

%% Collision with ground, where ground is r_y=0
for ts= 1:N
    v_x(ts+1) = v_x(ts) ;                    % Find new x velocity
    r_x(ts+1) = r_x(ts) + v_x(ts)*dt ;       % Find new x position
    r_y(ts+1) = r_y(ts) + v_y(ts)*dt ;       % Find new y position

    GPE(ts+1) = m*g*r_y(ts);                 % Gravitational Potential Energy calculated using m*g*h in J
    KE(ts+1)  = m*0.5*v_y(ts)^2;             % Kinetic Energy calculated using 0.5*m*v^2 in J
    TE(ts+1)  = GPE(ts)+KE(ts);              % Total Energy of Particle, sum of GPE and KE

    tt(ts+1)  = (ts+dt)*dt;
        if r_y(ts) <0                        % If condition that determines when projectile reaches the ground, where ground is y=0
            v_y(ts+1) =-Re* v_y(ts);         % Find new y velocity, as particle bounces back
            r_y(ts+1) = 0;                   % The  y position, as particle bounces back is 0
        else
            v_y(ts+1) = v_y(ts) + -g*dt;     % Find new y velocity, if no impact occurs
        end                                  % End If statement

             if TE(ts) < m*g                 % If condition determines if total energy of particle is less than 0.05
                 break                       % IF true, for loop is broken
             end                             % End If statement
end                                          % End for loop

%% Displaying Projectile Flight Information
hFigure = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.0, 0.0, 01, 01]);
set(gcf, 'Name', 'Projectile motion for a bouncing ball', 'NumberTitle', 'Off') 
subplot(3, 1, 1);
    plot(tt,r_y,'.')
    hold on 
    title("Projectile Motion: with particle Impact");
    legend("Euler's Method",'location','northeast')
    xlabel('Time(s)'); ylabel('Height(m)');
    go = gca; go.XAxisLocation = 'origin'; go.YAxisLocation = 'origin';
    hold off
subplot(3, 1, 2);
    plot(r_x,r_y,'.')
    hold on 
    legend("Euler's Method",'location','northeast')
    xlabel('Horizontal Distance(m)'); ylabel('Height(m)');
    go = gca; go.XAxisLocation = 'origin'; go.YAxisLocation = 'origin';
    hold off
subplot(3, 1, 3);
    hold on 
    plot(tt,GPE,'.','Color',[0.4660 0.6740 0.1880]);
    plot(tt,KE, '.','Color',[0.8500 0.3250 0.0980]);
    plot(tt,TE, '.-','Color',[0 0 0])
    title("Energy of Particle through out its motion");
    legend("GPE",'KE','TE','location','northeast')
    xlabel('Time(s)'); ylabel('Energy(J)');
    go = gca; go.XAxisLocation = 'origin'; go.YAxisLocation = 'origin';
    hold off

m1 = sprintf('\n      The total time of travel is %f seconds.\n', tt(:,end));
m2 = sprintf('      The total distance travelled by particle is %f meters.\n',r_x(:,end));
m3 = sprintf('      The maximum height reached by particle is %f meters.\n',max(r_y)); 
message = sprintf('%s', m1, m2,m3); disp(message)
%%
close all
    plot(tt,r_y,'.')
    hold on 
    title("Particle Impact: with dt=0.5");
    legend("Euler's Method",'location','northeast')
    xlabel('Time(s)'); ylabel('Height(m)');
    go = gca; go.XAxisLocation = 'origin'; go.YAxisLocation = 'origin';
    hold off