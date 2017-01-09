[samples, framerate] = audioread('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Sound Samples\Trumpet Violin\TA3_Vfs4.wav');

monoChannel = samples(:,1);

sampleLength = size(monoChannel);

windowSize = 2000;
hopSize = 1000;
nfft = 4000;

[S, F, T] = stft(monoChannel, windowSize, hopSize, nfft, framerate);
[x, t] = istft(S, hopSize, nfft, framerate);

audiowrite('C:\Users\Will\OneDrive\Uni\Individual Project\SoundSeparation\Output Sounds SMProto\istftTest.wav', x, framerate);