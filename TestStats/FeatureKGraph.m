k = [2, 3, 4, 5];
sinSDR = [-3.86812, 2.85519, -22.9901, -15.8334];
sinSIR = [41.738, 42.1597, 5.91709, 3.32598];
sinSAR = [-3.76851, 2.9057, -7.28391, -7.20458];

matSDR = [-5.37837, 2.49075, -6.21275, -7.23518];
matSIR = [2.29012, 18.7251, -4.9194, -4.18413];
matSAR = [-1.66888, 4.28889, 6.21312, 3.66141];

figure
subplot(3, 1, 1)
plot(k, sinSDR)
hold on
plot(k, matSDR)
ylim([-30, 10])
ylabel('SDR')
legend('Sinusoids', 'Matrix Factors')
title('Separation performance against number of sounds')

subplot(3, 1, 2)
plot(k, sinSIR)
hold on
plot(k, matSIR)
ylabel('SIR')

subplot(3, 1, 3)
plot(k, sinSAR)
hold on
plot(k, matSAR)
xlabel('Number of sounds')
ylabel('SAR')