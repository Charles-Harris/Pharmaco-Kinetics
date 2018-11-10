function hw4b

% Cleaning workplace
close all;
clear all;

% Recording randomized data
firstSmokerData = [0.401282	0.314772 0.236688 8.24E-02 5.09E-02 2.47E-02];
secondSmokerData = [0.316974 0.295283 0.230779 0.111714	8.17E-02 5.18E-02];
firstNonSmokerData = [0.333059 0.266863	0.268022 0.226398 0.251007 0.261069];
secondNonSmokerData = [0.37814 0.302747	0.307756 0.234712 0.185152 0.154593];

% Initializing values
D = 10; %mg
V = 25; %L

% Recorded time points
time = [1, 3, 5, 14, 18, 24]; % hours
timeInterval = (0:0.1:24); % hours

% Subject 1 - First Smoker

guess1(1) = .003; % First guess
guess1(2) = .004; % Second Guess
index = 2; 

while abs(newFunction(guess1(index),time,firstSmokerData)) > .01
    guess1(index+1) = guess1(index) - newFunction(guess1(index),time,firstSmokerData)*(guess1(index)-guess1(index-1))/((newFunction(guess1(index), time, firstSmokerData)-newFunction(guess1(index-1), time, firstSmokerData)));
    index = index + 1;
end

kValueFirstSmoker = guess1(index);
firstSolution = (D/V) * exp(-kValueFirstSmoker*timeInterval);

% Subject 2 - Second Smoker

guess2(1) = .01; % First guess
guess2(2) = .02; % Second guess
index = 2;

while abs(newFunction(guess2(index), time, secondSmokerData))> .01
    guess2(index+1) = guess2(index) - (newFunction(guess2(index),time,secondSmokerData))*(guess2(index)-guess2(index-1))/(newFunction(guess2(index),time,secondSmokerData)-newFunction(guess2(index-1),time,secondSmokerData));
    index = index+1;
end

kValueSecondSmoker = guess2(index);
secondSolution = (D/V)*exp(-kValueSecondSmoker*timeInterval);

% Subject 3 - First NonSmoker

guess3(1) = .05; % First guess
guess3(2) = .06; % Second guess
index = 2;

while abs(newFunction(guess3(index),time,firstNonSmokerData)) > .01
    guess3(index+1) = guess3(index) - newFunction(guess3(index),time,firstNonSmokerData)*(guess3(index)-guess3(index-1))/(newFunction(guess3(index),time,firstNonSmokerData) - newFunction(guess3(index-1),time,firstNonSmokerData));
    index = index + 1;
end

kValueFirstNonSmoker = guess3(index);
thirdSolution=(D/V)*exp(-kValueFirstNonSmoker*timeInterval);

% Subject 4 - Second Non Smoker

guess4(1) = .03; % First guess
guess4(2) = .04; % Second guess
index = 2;

while abs(newFunction(guess4(index),time,secondNonSmokerData)) > .01
    guess4(index+1) = guess4(index)-newFunction(guess4(index),time,secondNonSmokerData)*(guess4(index)-guess4(index-1))/(newFunction(guess4(index),time,secondNonSmokerData)-newFunction(guess4(index-1),time,secondNonSmokerData));
    index = index+1;
end

kValueSecondNonSmoker = guess4(index);
fourthSolution = (D/V)*exp(-kValueSecondNonSmoker*timeInterval)

% Making the plots

subplot(2,2,1);
plot(time, firstSmokerData,'ob');
hold on;
plotOne = plot(timeInterval,firstSolution,'b')
title('Smoker 1 Concentration Values');
xlabel('Time (h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legend1=strcat('K Est:', num2str(kValueFirstSmoker));
legend(plotOne,legend1);

subplot(2,2,2);
plot(time, secondSmokerData,'or');
hold on;
plotTwo = plot(timeInterval,secondSolution,'r');
title('Smoker 2 Concentration Values');
xlabel('Time (h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legend2=strcat('K Est:', num2str(kValueSecondSmoker));
legend(plotTwo,legend2);

subplot(2,2,3);
plot(time, firstSmokerData,'oy');
hold on;
plotThree = plot(timeInterval,thirdSolution,'y');
title('NonSmoker 1 Concentration Values');
xlabel('Time (h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legend3=strcat('K Est:', num2str(kValueFirstNonSmoker));
legend(plotThree,legend3);

subplot(2,2,4);
plot(time, secondNonSmokerData,'og');
hold on;
plotFour = plot(timeInterval,fourthSolution,'g')
title('NonSmoker 2 Concentration Values');
xlabel('Time (h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legend4=strcat('K Est:', num2str(kValueFirstSmoker));
legend(plotFour,legend4);

% Secant method
function f = newFunction(x,time,data);
Y = data.*time.*exp(time.*-x);
firstSum = sum(Y);
Z = time.*exp(time.*-2*x);
secondSum = sum(Z);
f = (firstSum-(D/V)*secondSum);

end

end









