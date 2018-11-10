% This program will solve a series of differential equations using Matlab's
% built-in ODE solver. 
% This program will solve the problem presented in part 2, but will also
% account for additional concentrations of nicotine accumulated from
% smoking cigarettes

function hw2c

clear all;
close all;

% iii. First Nicotine Pulse (From time t=0 to time t=1 (in minutes))

% Initial Concentrations
initialL = 11; % micromoles
initialR = 11; % micromoles
initialD = 11; % micromoles
initialLR = 0; % micromoles
initialLD = 0; % micromoles

% Setup and Call to ODE45
trange0 = [0 1];  % time interval 0 in minutes
x0=[initialL+15 initialR initialD initialLR initialLD]; % initial condition
[tpulse0,xpulse0]=ode45(@ODE45SPde,trange0, x0,odeset('RelTol', 1e-5, 'AbsTol', 1e-5));


% iv. Second Nicotine Pulse (From time t=1 to time t=2 (in minutes))

trange1 = [1 2]; %time interval 1 in minutes
x1 =[xpulse0(end,1)+15 xpulse0(end,2) xpulse0(end,3) xpulse0(end,4) xpulse0(end,5)];
[tpulse1,xpulse1]=ode45(@ODE45SPde, trange1, x1, odeset('RelTol',1e-5,'AbsTol',1e-5));


% v. Third Nicotine Pulse (From time t=2 to time t=3 (in minutes))

trange2 = [2 3];
x2=[xpulse1(end,1)+15 xpulse1(end,2) xpulse1(end,3) xpulse1(end,4) xpulse1(end,5)];
[tpulse2,xpulse2]=ode45(@ODE45SPde, trange2, x2, odeset('RelTol',1e-5,'AbsTol',1e-5));

% vi. Vector concatenations

timeConcatentaion = [tpulse0; tpulse1; tpulse2];
concentrationConcatentation = [xpulse0; xpulse1; xpulse2];

%
%Plotting

%Subplot 1
subplot(2,2,1)
plot(timeConcatentaion,concentrationConcatentation(:,2),'r-')
xlim([-.5,3])
ylim([0, 13])
xlabel('Time (min)')
ylabel('Concentration of R (uM)')
title('Nicotine Pulse - [R] (uM) vs Time (min)')

%Subplot 2
subplot(2,2,2)
plot(timeConcatentaion,concentrationConcatentation(:,4), 'b-')
xlim([-.5,3])
xlabel('Time (min)')
ylabel('Concentration of LR (uM)')
title('Nicotine Pulse - [LR] (uM) vs Time (min)')

%Subplot 3
subplot(2,2,3)
plot(timeConcatentaion,concentrationConcatentation(:,3),'g-')
xlim([-.5,3])
ylim([0,13])
xlabel('Time (min)')
ylabel('Concentration of D (uM)')
title('Nicotine Pulse - [D] (uM) vs Time (min)')

%Subplot 4
subplot(2, 2, 4)
plot(timeConcatentaion,concentrationConcatentation(:,5), 'm-')
xlim([-.5,3])
xlabel('Time (min)')
ylabel('Concentration of LD (uM)')
title('Nicotine Pulse - [LD] (uM) vs Time (min)')

return

% ODE45SPde function:
function dxdt= ODE45SPde(t,x)

% Values of Rate constants
k1 = 600; % uM^-1min^-1
k2 = 6000; % min^-1
k3 = 28.8; % min^-1
k4 = 0.288; % min^-1
k5 = 600; % uM^-1min^-1
k6 = 60; % min^-1
k7 = 0.8; % uM^-1min^-1
k8 = 0.8; % uM^-1min^-1
k9 = 36.0; % uMmin^-1
k10 = 120.0; %min^-1

dxdt=[-k1*x(1)*x(2)+k2*x(4)-k5*x(1)*x(3)+k6*x(5)+k9-k10*x(1);
      -k1*x(1)*x(2)+k2*x(4)-k7*x(1)*x(2)+k8*x(1)*x(3);
      -k5*x(1)*x(3)+k6*x(5)+k7*x(1)*x(2)-k8*x(1)*x(3);
       k1*x(1)*x(2)-k2*x(4)-k3*x(4)+k4*x(5);
       k3*x(4)-k4*x(5)+k5*x(1)*x(3)-k6*x(5)];
return