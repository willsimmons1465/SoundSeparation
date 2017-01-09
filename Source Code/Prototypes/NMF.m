WINDOWSIZE = 1000;
HOPSIZE = 500;
NFFT = 20000;
NUMSOUNDS = 2;
FILE = 'Guitar Trumpet\GD3vlfn_TGs41fn';

%Get input audio
inputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Sound Samples\\%s.wav', FILE);
[samples, framerate] = audioread(inputFilename);

clear inputFilename;

%Perform STFT
[spectrum, ~, ~] = stft(samples, WINDOWSIZE, HOPSIZE, NFFT, framerate);

clear samples;

%Perform NMF
absSpectrum = abs(spectrum);
[mix, source] = nnmf(absSpectrum,NUMSOUNDS);

%Generate and output sounds
for i = 1:NUMSOUNDS
    soundSpectrum = mix(:,i) * source(i,:);
    [soundSamples, ~] = istft(soundSpectrum, HOPSIZE, NFFT, framerate);
    outputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Output Sounds NMFProto\\%s_%d.wav', FILE, i);
    audiowrite(outputFilename, soundSamples, framerate);
end