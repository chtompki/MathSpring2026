daysAfterBirth  = [0 6 10 13 17 20 28];
weightYoungLeaves = [7 18 45 40 32 30.5 30];
weightMatureLeaves = [7 16 19 15 12 10.5 10];

youngWeightNatural = naturalCubicSplines(daysAfterBirth,weightYoungLeaves,3)
matureWeightNatural = naturalCubicSplines(daysAfterBirth,weightMatureLeaves,3)
meanWeightNatural = (youngWeightNatural + matureWeightNatural)/2

youngWeightMatlab = spline(daysAfterBirth,weightYoungLeaves,3)
matureWeightMatlab = spline(daysAfterBirth,weightMatureLeaves,3)
meanWeightMatlab = (youngWeightMatlab + matureWeightMatlab)/2

% >> gl3prob6
%
% youngWeightNatural =
%
%     7.1888
%
%
% matureWeightNatural =
%
%    11.5619
%
%
% meanWeightNatural =
%
%     9.3754
%
%
% youngWeightMatlab =
%
%    -2.7094
%
%
% matureWeightMatlab =
%
%    10.4508
%
%
% meanWeightMatlab =
%
%     3.8707

xstar=0:0.25:30;
plottableYoungWeightNatural = naturalCubicSplines(daysAfterBirth,weightYoungLeaves,xstar);
plottableMatureWeightNatural = naturalCubicSplines(daysAfterBirth,weightMatureLeaves,xstar);
meanWeights = (weightYoungLeaves + weightMatureLeaves)./ 2;
plottableMeanWeightsNatural = naturalCubicSplines(daysAfterBirth,meanWeights,xstar);

figure
hold on
plot(xstar,plottableYoungWeightNatural,'k-');
plot(xstar, plottableMatureWeightNatural, 'r-');
plot(xstar, plottableMeanWeightsNatural, 'k--')
xlabel('Days After Birth');
ylabel('Weight');
legend('Young Leaves', 'Mature Leaves', 'Mean Weights');
hold off;

