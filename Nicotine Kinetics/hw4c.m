function hw4c

% Recording randomized data
firstSmokerData = [0.401282	0.314772 0.236688 8.24E-02 5.09E-02 2.47E-02];
secondSmokerData = [0.316974 0.295283 0.230779 0.111714	8.17E-02 5.18E-02];
firstNonSmokerData = [0.333059 0.266863	0.268022 0.226398 0.251007 0.261069];
secondNonSmokerData = [0.37814 0.302747	0.307756 0.234712 0.185152 0.154593];

% Recorded time points
time = [1, 3, 5, 14, 18, 24]; % hours
timeInterval = (0:0.1:24); % hours

% Initial Dose
D = 10; %mg

% Vector of guesses
firstKValuesGuess = [.3 .3];

% Nelder-Mead Method
[kValueSmoker1, result] = fminsearch(@(KValues) fmin(KValues, D, time, firstSmokerData), firstKValuesGuess);
[kValueSmoker2, result] = fminsearch(@(KValues) fmin(KValues, D, time, secondSmokerData), firstKValuesGuess);
[kValueNonSmoker1, result] = fminsearch(@(KValues) fmin(KValues, D, time, firstNonSmokerData), firstKValuesGuess);
[kValueNonSmoker2, result] = fminsearch(@(KValues) fmin(KValues, D, time, secondNonSmokerData), firstKValuesGuess);

% Putting the K Values into a vector
KValues = vertcat(kValueSmoker1, kValueSmoker2, kValueNonSmoker1,kValueNonSmoker2);

% Assigning K Values
smokerOneK = KValues(1,1);
smokerTwoK = KValues(2,1);
nonSmokerOneK = KValues(3,1);
nonSmokerTwoK = KValues(4,1);

% Assigning V Values
smokerOneV = KValues(1,2);
smokerTwoV = KValues(2,2);
nonSmokerOneV = KValues(3,2);
nonSmokerTwoV = KValues(4,2);

% Setting up the graphs

% Plot One
Y1 = (D/smokerOneV)*exp(-smokerOneK*timeInterval);
subplot(2,2,1)
plotOne = plot(timeInterval, Y1, 'b');
hold on;
plot(time, firstSmokerData,'ob');
hold on;
title('Smoker 1 Concentration v Time');
xlabel('Time(h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legendOne = strcat('K: ', num2str(smokerOneK), ' V: ', num2str(smokerOneV));
legend(plotOne,legendOne);

% Plot Two
Y2 = (D/smokerTwoV)*exp(-smokerTwoK*timeInterval);
subplot(2,2,2)
plotTwo = plot(timeInterval, Y2, 'r');
hold on;
plot(time, secondSmokerData,'or');
hold on;
title('Smoker 2 Concentration v Time');
xlabel('Time(h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legendTwo = strcat('K: ', num2str(smokerTwoK), ' V: ', num2str(smokerTwoV));
legend(plotTwo,legendTwo);

% Plot Three
Y3 = (D/nonSmokerOneV)*exp(-nonSmokerOneK*timeInterval);
subplot(2,2,3)
plotThree = plot(timeInterval, Y3, 'y');
hold on;
plot(time, firstNonSmokerData,'oy');
hold on;
title('NonSmoker 1 Concenttration v Time');
xlabel('Time(h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legendThree = strcat('K: ', num2str(nonSmokerOneK), ' V: ', num2str(nonSmokerOneV));
legend(plotThree,legendThree);

% Plot Four
Y4 = (D/nonSmokerTwoV)*exp(-nonSmokerTwoK*timeInterval);
subplot(2,2,4)
plotFour = plot(timeInterval, Y4, 'm');
hold on;
plot(time, secondNonSmokerData,'om');
hold on;
title('NonSmoker 2 Concentration v Time');
xlabel('Time(h)');
ylabel('Concentration (mg/L)');
axis([0 24 0 0.8]);
legendFour = strcat('K: ', num2str(nonSmokerTwoK), ' V: ', num2str(nonSmokerTwoV));
legend(plotFour,legendFour);

% Sum of Squares Function
function f = fmin(KValues, D, time, data)
K = KValues(1);
V = KValues(2);
Z = ((data - (D/V).*exp(time.*-K)).^2);
sumOne = sum(Z);
f = sumOne;
end

end




