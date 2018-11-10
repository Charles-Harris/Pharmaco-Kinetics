% The objective of this program is to analyze the receptor ligand kinetics of the
% ?2?4 nicotinic acetylcholine receptor
% This program will approximate a series of differential equations using
% Euler's method

function hw2a

clear all;
close all;


initialTime = 0; %minutes
finalTime = 1/60; % minutes
h=0.000001; %Integration step size % minutes
t=[initialTime:h:finalTime]; %time vector
n=length(t);

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

% Storing Initial Conditions
x(1,:)=[initialL initialR initialD initialLR initialLD]; %Initial conditions stored in first row of x

%Propagation loop for Euler's method. Solution stored in subsequent rows
for i=1:n-1
    x(i+1,:)=x(i,:)+h*SPde(x(i,:),t(i))';
end

%Plotting
subplot(2,2,1)
plot(t,x(:,2),'r-')
xlabel('Time (min)')
ylabel('Concentration of R (uM)')
title('Concentration of R (uM) vs Time (min)')
subplot(2,2,2)
plot(t,x(:,4), 'b-')
xlabel('Time (min)')
ylabel('Concentration of LR (uM)')
title('Concentration of LR (uM) vs Time (min)')
subplot(2,2,3)
plot(t,x(:,3),'g-')
xlabel('Time (min)')
ylabel('Concentration of D (uM)')
title('Concentration of D (uM) vs Time (min)')
subplot(2, 2, 4)
plot(t,x(:,5), 'm-')
xlabel('Time (min)')
ylabel('Concentration of LD (uM)')
title('Concentration of LD (uM) vs Time (min)')

return

% SPde function:
function dxdt= SPde(x,t)

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