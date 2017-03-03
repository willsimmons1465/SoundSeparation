WINDOWSIZE = 1000;
HOPSIZE = 500;
NFFT = 2000;
%FILE = 'OptimiseTests\saxophone_E5_1_forte_normal';
%FILE = 'OptimiseTests\trumpet_Gs4_1_forte_normal';
%FILE = 'OptimiseTests\SE51fn_TGs41fn_0';
FILE = 'OptimiseTests\SE51fn_TGs41fn_1';

%inputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Sound Samples\\%s.wav', FILE);
inputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Output Sounds SM\\%s.wav', FILE);
[samples, framerate] = audioread(inputFilename);

if size(samples, 2) == 1
    samples = [samples, samples];
end

paddedSamples = zeros(104832, 2);
paddedSamples(1:size(samples, 1), :) = samples;

[Lspectrum, f, t] = stft(paddedSamples(:,1), WINDOWSIZE, HOPSIZE, NFFT, framerate);
[Rspectrum, ~, ~] = stft(paddedSamples(:,2), WINDOWSIZE, HOPSIZE, NFFT, framerate);

spectrum = sqrt(real(Lspectrum).^2 + imag(Lspectrum).^2 + real(Rspectrum).^2 + imag(Rspectrum).^2);
logSpec = log2(spectrum);
imagesc([min(t), max(t)], [min(f), max(f)], logSpec, [-6, max(max(logSpec))]);
%title('Spectrogram of original saxophone sample');
%title('Spectrogram of original trumpet sample');
%title('Spectrogram of reconstructed saxophone sample');
title('Spectrogram of reconstruced trumpet sample');
xlabel('Time/s');
ylabel('Frequency/Hz');