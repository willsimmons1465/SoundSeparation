inSNR = [-14.1935, -14.1805, -13.2653, -12.5096, -11.6817, -8.1599, -2.1393, 3.8813, 11.8401];
outSNR = [-365.063, -352.602, -81.4467, 9.31435, 35.9491, 45.9647, 54.0935, 55.1898, 57.9213];
outSNREr = [477.789, 469.126, 127.597, 55.3973, 18.9779, 8.21627, 7.59442, 6.89302, 8.4362];
inSNRMat = [-14.1805, -11.6817, -8.1599, -2.1393, 3.8813, 11.8401];
matSNR = [-14.9844, -11.8032, -7.03229, 1.58634, 9.24763, 17.7985];
matSNREr = [5.03715, 4.9634, 5.12386, 5.76312, 6.55015, 6.77564];

%errorbar(inSNR, outSNR, outSNREr)
plot(inSNR, outSNR)
hold on
plot(inSNRMat, matSNR)
xlabel('Mean SNR of inputs/dB')
ylabel('Mean SNR of outputs/dB')
legend('Sin', 'NMF', 'Location', 'southeast')
%title('Noise reduction')