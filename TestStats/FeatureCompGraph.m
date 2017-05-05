featureCounts = [10, 12, 19, 22, 24, 28, 47, 72, 130, 179, 221];
featureCounts = featureCounts ./ 2;
sinSDR = [-17.5519, -15.2417, -6.63394, -5.33899, -5.11012, -5.03388, -8.23324, -10.3056, -12.0895, -8.48039, -7.52456];
sinSDREr = [20.9057, 19.9599, 12.3659, 11.0721, 10.9822, 11.1301, 15.2445, 17.63, 14.5556, 11.8642, 13.917];
sinSIR = [21.4367, 24.5176, 34.9656, 37.5821, 36.0798, 33.577, 31.2791, 29.6929, 22.0113, 26.7663, 28.5192];
sinSIREr = [30.9733, 30.708, 26.2405, 25.0639, 25.6171, 23.9871, 28.0367, 26.7259, 20.885, 19.0216, 17.3623];
sinSAR = [-1.04206, -0.69223, -0.155166, -0.297683, -0.102584, -0.0484372, -1.21733, -1.29858, -3.78852, -4.13248, -3.38612];
sinSAREr = [6.04595, 5.49339, 5.13122, 5.699, 5.62412, 5.27561, 4.65396, 6.41558, 6.59486, 6.47318, 8.37485];

matSDR = [1.27617, 0.349679, 0.653252, 0.955421, -0.43112, 0.33252, -0.618318, 0.00160657, 0.357906, -0.23749, 0.13676];
matSDREr = [5.42448, 5.44717, 5.33583, 5.63631, 5.21494, 4.76763, 5.13074, 4.18638, 4.27786, 4.29291, 4.5345];
matSIR = [10.0806, 7.07855, 8.15148, 9.33958, 6.10001, 5.85867, 5.7812, 4.79253, 5.91472, 4.11621, 6.06095];
matSIREr = [9.40915, 7.52763, 8.87433, 9.52119, 8.72099, 6.14115, 6.84277, 6.52145, 6.46788, 6.35826, 8.48208];
matSAR = [6.28145, 6.56551, 7.01717, 6.66738, 7.4615, 8.29515, 7.13291, 8.34871, 8.08928, 8.17739, 8.05827];
matSAREr = [3.57425, 3.556, 3.53345, 4.06359, 4.41694, 4.1844, 3.50765, 3.6195, 3.46184, 3.21701, 3.56521];

figure
subplot(3, 1, 1)
%errorbar(featureCounts, sinSDR, sinSDREr)
plot(featureCounts, sinSDR)
hold on
%errorbar(featureCounts, matSDR, matSDREr)
plot(featureCounts, matSDR)
ylabel('SDR/dB')
%legend('Sinusoids', 'Matrix Factors')
legend('Sin', 'NMF')
%title('Separation performance against number of features extracted')

subplot(3, 1, 2)
%errorbar(featureCounts, sinSIR, sinSIREr)
plot(featureCounts, sinSIR)
hold on
%errorbar(featureCounts, matSIR, matSIREr)
plot(featureCounts, matSIR)
ylabel('SIR/dB')

subplot(3, 1, 3)
%errorbar(featureCounts, sinSAR, sinSAREr)
plot(featureCounts, sinSAR)
hold on
%errorbar(featureCounts, matSAR, matSAREr)
plot(featureCounts, matSAR)
xlabel('Number of features')
ylabel('SAR/dB')