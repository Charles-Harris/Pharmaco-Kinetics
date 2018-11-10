% This program will solve a series of differential equations using Matlab's
% built-in ODE solver. 
% The equations being solved are the same ones from Part 2

function hw2b

clear all;
close all;

initialTime = 0; %minutes
finalTime = 1/60; % minutes

% Initial Concentrations
initialL = 11; % micromoles
initialR = 11; % micromoles
initialD = 11; % micromoles
initialLR = 0; % micromoles
initialLD = 0; % micromoles

% Values of Rate constants
% k1 = 600; %uM^-1min^-1
% k2 = 6000; % min^-1
% k3 = 28.8; % min^-1
% k4 = 0.288; % min^-1
% k5 = 600; %uM^-1min^-1
% k6 = 60; % min^-1
% k7 = 0.8; %uM^-1min^-1
% k8 = 0.8; %uM^-1min^-1

% Commented out here because they arent used until the next function

% Setup and Call to ODE45
trange =[initialTime finalTime];  % time interval
x0=[initialL initialR initialD initialLR initialLD]; % initial condition
[t,x]=ode45(@ODE45SPde,trange, x0,odeset('RelTol', 1e-5, 'AbsTol', 1e-5));

%Plotting
subplot(2,2,1)
plot(t,x(:,2),'r-')
xlabel('Time (min)')
ylabel('Concentration of R (uM)')
title('ODE45 - [R] (uM) vs Time (min)')
subplot(2,2,2)
plot(t,x(:,4), 'b-')
xlabel('Time (min)')
ylabel('Concentration of LR (uM)')
title('ODE45 - [LR] (uM) vs Time (min)')
subplot(2,2,3)
plot(t,x(:,3),'g-')
xlabel('Time (min)')
ylabel('Concentration of D (uM)')
title('ODE45 - [D] (uM) vs Time (min)')
subplot(2, 2, 4)
plot(t,x(:,5), 'm-')
xlabel('Time (min)')
ylabel('Concentration of LD (uM)')
title('ODE45 - [LD] (uM) vs Time (min)')

return

% ODE45SPde function:
function dxdt= ODE45SPde(t,x)

% Values of Rate constants
k1 = 600; %uM^-1min^-1
k2 = 6000; % min^-1
k3 = 28.8; % min^-1
k4 = 0.288; % min^-1
k5 = 600; %uM^-1min^-1
k6 = 60; % min^-1
k7 = 0.8; %uM^-1min^-1
k8 = 0.8; %uM^-1min^-1

dxdt=[-k1*x(1)*x(2)+k2*x(4)-k5*x(1)*x(3)+k6*x(5);
      -k1*x(1)*x(2)+k2*x(4)-k7*x(1)*x(2)+k8*x(1)*x(3);
      -k5*x(1)*x(3)+k6*x(5)+k7*x(1)*x(2)-k8*x(1)*x(3);
       k1*x(1)*x(2)-k2*x(4)-k3*x(4)+k4*x(5);
       k3*x(4)-k4*x(5)+k5*x(1)*x(3)-k6*x(5)];
return