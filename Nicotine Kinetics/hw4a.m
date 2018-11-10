function hw4a

% cleaning workspace
clear all
close all

% Initializing variables

% Initializing time range
trange = [0,8];

% Initializing concentrations
x0 = [50 0 0];
x1 = [100 0 0];

% Setting up ODE45
%%%%%%%%%%%%%%%%%%%%%%

% Males 50 g Dose
[t1,xm1] = ode45(@ODE45SPde1, trange , x0, odeset('RelTol',1e-5, 'AbsTol', 1e-5));

%Males 100 g Dose
[t2,xm2] = ode45(@ODE45SPde1, trange , x1, odeset('RelTol',1e-5, 'AbsTol', 1e-5));

%Females 50 g Dose
[t3,xf1] = ode45(@ODE45SPde2, trange , x0, odeset('RelTol',1e-5, 'AbsTol', 1e-5));

%Females 100 g Dose
[t4,xf2] = ode45(@ODE45SPde2, trange , x1, odeset('RelTol',1e-5, 'AbsTol', 1e-5));

% Plotting
subplot(2,1,1);
plot(t1,xm1(:,1),'b');
hold on;
plot(t2, xm2(:,1),'g')
hold on;
plot(t3,xf1(:,1),'y');
hold on;
plot(t4,xf2(:,1),'m');
xlim([0,1])
legend('Males (50g)', 'Males (100g)', 'Females (50g)', 'Females (100g)');
xlabel('Time (h)')
ylabel('Mass of Alcohol Present (g)')
title ('Graph of Mass of Alcohol Present in the Stomach (g) vs Time (h)')
subplot(2,1,2);
plot(t1,xm1(:,3)/420, 'b');
hold on;
plot(t2,xm2(:,3)/420,'g');
hold on;
plot(t3,xf1(:,3)/300,'y');
hold on;
plot(t4,xf2(:,3)/300,'m');
legend('Males (50g)', 'Males (100g)', 'Females(50g)', 'Females(100g)');
xlabel('Time(h)');
ylabel('Blood Alcohol Content (% body mass)')
title('Graph of BAC (% mass) vs Time (h)');
figure;

plot(t1,xm1(:,3)/420, 'b');
hold on;
plot(t2,xm2(:,3)/420,'g');
hold on;
plot(t3,xf1(:,3)/300,'y');
hold on;
plot(t4,xf2(:,3)/300,'m');
hold on;
plot(t1,.08,'k');
legend('Males (50g)', 'Males (100g)', 'Females(50g)', 'Females(100g)','Legal Driving Concentration (.08% BAC)');
xlabel('Time(h)');
ylabel('Blood Alcohol Content (% body mass)')
title('Graph of BAC (% mass) vs Time (h)');

return

function dxdt1 = ODE45SPde1(t,x)

%males
ks = 10; % h^-1
ka = 11000; % h^-1
Vm = 45; % g/h
Km = 50; % g

dxdt1 = [-ks*x(1);
        ks*x(1)-ka*x(2);
        ka*x(2)-((Vm*x(3))/(Km+x(3)))];
    
 return
 
 function dxdt2 = ODE45SPde2(t,x)

% females
ks = 10; % h^-1
ka = 11000; % h^-1
Vm = 40; % g/h
Km = 50; % g

dxdt2 = [-ks*x(1);
        ks*x(1)-ka*x(2);
        ka*x(2)-((Vm*x(3))/(Km+x(3)))];
    
 return
 



