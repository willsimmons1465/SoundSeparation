featureCounts = [28, 47, 72, 130, 179, 221];
featureCounts = featureCounts ./ 2;
sinSDR = [-5.03388, -8.23324, -10.3056, -12.0895, -8.48039, -7.52456];
sinSDREr = [11.1301, 15.2445, 17.63, 14.5556, 11.8642, 13.917];
sinSIR = [33.577, 31.2791, 29.6929, 22.0113, 26.7663, 28.5192];
sinSIREr = [23.9871, 28.0367, 26.7259, 20.885, 19.0216, 17.3623];
sinSAR = [-0.0484372, -1.21733, -1.29858, -3.78852, -4.13248, -3.38612];
sinSAREr = [5.27561, 4.65396, 6.41558, 6.59486, 6.47318, 8.37485];

matSDR = [0.33252, -0.618318, 0.00160657, 0.357906, -0.23749, 0.13676];
matSDREr = [4.76763, 5.13074, 4.18638, 4.27786, 4.29291, 4.5345];
matSIR = [5.85867, 5.7812, 4.79253, 5.91472, 4.11621, 6.06095];
matSIREr = [6.14115, 6.84277, 6.52145, 6.46788, 6.35826, 8.48208];
matSAR = [8.29515, 7.13291, 8.34871, 8.08928, 8.17739, 8.05827];
matSAREr = [4.1844, 3.50765, 3.6195, 3.46184, 3.21701, 3.56521];

figure
subplot(3, 1, 1)
errorbar(featureCounts, sinSDR, sinSDREr)
%plot(featureCounts, sinSDR)
hold on
errorbar(featureCounts, matSDR, matSDREr)
%plot(featureCounts, matSDR)
ylabel('SDR')
legend('Sinusoids', 'Matrix Factors')
title('Separation performance against number of features extracted')

subplot(3, 1, 2)
errorbar(featureCounts, sinSIR, sinSIREr)
%plot(featureCounts, sinSIR)
hold on
errorbar(featureCounts, matSIR, matSIREr)
%plot(featureCounts, matSIR)
ylabel('SIR')

subplot(3, 1, 3)
errorbar(featureCounts, sinSAR, sinSAREr)
%plot(featureCounts, sinSAR)
hold on
errorbar(featureCounts, matSAR, matSAREr)
%plot(featureCounts, matSAR)
xlabel('Number of features')
ylabel('SAR')