inSNR = [-14.1805, -11.6817, -8.1599, -2.1393, 3.8813, 11.8401];
outSNR = [-352.602, 35.9491, 45.9647, 54.0935, 55.1898, 57.9213];
outSNREr = [469.126, 18.9779, 8.21627, 7.59442, 6.89302, 8.4362];
matSNR = [-14.9844, -11.8032, -7.03229, 1.58634, 9.24763, 17.7985];
matSNREr = [5.03715, 4.9634, 5.12386, 5.76312, 6.55015, 6.77564];

%errorbar(inSNR, outSNR, outSNREr)
plot(inSNR, outSNR)
hold on
plot(inSNR, matSNR)
xlabel('Mean SNR of input')
ylabel('Mean SNR of outputs')
legend('Sinusoids', 'Matrix Factors', 'Location', 'southeast')
title('Noise reduction')