WINDOWSIZE = 1000;
HOPSIZE = 500;
NFFT = 2000;
%THRESHOLD = 8;
THRESHOLD = 0.01;
TRAJECTORYDELTAFMAX = 0.02;
NUMSOUNDS = 2;
FILE = 'Guitar Trumpet\GD3vlfn_TGs41fn';

rng default;


%Get input audio
inputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Sound Samples\\%s.wav', FILE);
[samples, framerate] = audioread(inputFilename);

clear inputFilename;

%Force stereo handling
if size(samples, 2) == 1
    samples = [samples, samples];
end

%Perform STFT
%[Lspectrum, ~, ~] = stft(samples(:,1), WINDOWSIZE, HOPSIZE, NFFT, framerate);
%[Rspectrum, ~, ~] = stft(samples(:,2), WINDOWSIZE, HOPSIZE, NFFT, framerate);
wname = 'bump';
Lspectrum = cwt(samples(:,1),wname,framerate);
Rspectrum = cwt(samples(:,2),wname,framerate);

clear samples;

%Threshold
thresholdedSpectrum = sqrt(real(Lspectrum).^2 + imag(Lspectrum).^2 + real(Rspectrum).^2 + imag(Rspectrum).^2);
thresholdedSpectrum(thresholdedSpectrum < THRESHOLD) = 0;

imagesc(thresholdedSpectrum);

%Extract sinusoids
[peaks, peakLocations] = findpeaks(thresholdedSpectrum(:,1));
sinAmps = peaks(:);
sinLSpec = Lspectrum(peakLocations, 1);
sinRSpec = Rspectrum(peakLocations, 1);
sinFreqBins = peakLocations(:);
for i = 2:size(thresholdedSpectrum,2)
    [peaks, peakLocations] = findpeaks(thresholdedSpectrum(:,i));
    numPeaks = size(peaks(:));
    numSins = size(sinFreqBins(:,i-1));
    if numPeaks(1) == 0
        sinAmps(:, i) = 0;
        sinLSpec(:, i) = 0;
        sinRSpec(:, i) = 0;
        sinFreqBins(:, i) = 0;
    end
    
    for j = 1:numPeaks(1)
        peakLoc = peakLocations(j);
        peak = peaks(j);
        matches = find(abs(log(sinFreqBins(:,i-1)/peakLoc)) < log(1+TRAJECTORYDELTAFMAX));
        if isempty(matches)
            sinAmps(numSins(1)+1,i) = peak;
            sinLSpec(numSins(1)+1,i) = Lspectrum(peakLoc,i);
            sinRSpec(numSins(1)+1,i) = Rspectrum(peakLoc,i);
            sinFreqBins(numSins(1)+1,i) = peakLoc;
            numSins = size(sinFreqBins(:,i-1));
        else
            sinAmps(matches, i) = peak;
            sinLSpec(matches, i) = Lspectrum(peakLoc,i);
            sinRSpec(matches, i) = Rspectrum(peakLoc, i);
            sinFreqBins(matches, i) = peakLoc;
        end
    end
end

clear thresholdedSpectrum i j matches numPeaks numSins peak peakLoc peakLocations peaks;

%Normalise sinusoid data
meanAmps = MeanIgnoringZeros(sinAmps);
meanFreqs = MeanIgnoringZeros(sinFreqBins);
normAmps = bsxfun(@rdivide, sinAmps, meanAmps);
normFreqs = bsxfun(@rdivide, sinFreqBins, meanFreqs);

%Identify stereo positions
sinStereo = (sinLSpec.^2) ./ ((sinLSpec.^2) + (sinRSpec.^2));
sinStereo(isnan(sinStereo)) = 0.5;

%START Distance Matrix
numSins = size(sinAmps, 1);
MinFreq = min(meanFreqs);

df = zeros(numSins, numSins);
da = zeros(numSins, numSins);
dh = Inf(numSins, numSins);

for i = 1:numSins
    for j = i:numSins
        df(i,j) = sum((normFreqs(i,:) - normFreqs(j,:)).^2);
        da(i,j) = sum((normAmps(i,:) - normAmps(j,:)).^2);
        
        maxB = floor(meanFreqs(j)*(1+TRAJECTORYDELTAFMAX)/MinFreq);
        
        for a = 1:floor(meanFreqs(i)*(1+TRAJECTORYDELTAFMAX)/MinFreq)
            b = a*meanFreqs(j)/meanFreqs(i);
            if b >= maxB+1
                break;
            end
            
            if b >= 1
                dhab = abs(log((meanFreqs(i)*floor(b))/(meanFreqs(j)*a)));
                if dhab < dh(i,j)
                    dh(i,j) = dhab;
                end
            end
            
            if b > maxB
                break;
            end
            
            dhab = abs(log((meanFreqs(i)*ceil(b))/(meanFreqs(j)*a)));
            if dhab < dh(i,j)
                dh(i,j) = dhab;
            end
        end
        
        df(j,i) = df(i,j);
        da(j,i) = da(i,j);
        dh(j,i) = dh(i,j);
    end
end

WF = 1;
WA = 2.5;
WH = 1000;
d = WF*df + WA*da + WH*dh;
%END Distance Matrix

%Cluster
%featureVectors = [normAmps normFreqs sinStereo];
featureVectors = [d sinStereo];
clusterIndices = kmeans(featureVectors, NUMSOUNDS);

clear featureVectors sinAmps sinLSpec sinRSpec sinStereo;

for i = 1:NUMSOUNDS
    soundFreqBins = sinFreqBins(clusterIndices == i, :);
    soundLSpec = zeros(size(Lspectrum));
    soundRSpec = zeros(size(Rspectrum));
    for t = 1:size(Lspectrum,2)
        nonzeroFreqs = soundFreqBins(soundFreqBins(:,t) ~= 0, t);
        soundLSpec(nonzeroFreqs,t) = Lspectrum(nonzeroFreqs,t);
        soundRSpec(nonzeroFreqs,t) = Rspectrum(nonzeroFreqs,t);
    end
    %[soundLsamples, ~] = istft(soundLSpec, HOPSIZE, NFFT, framerate);
    %[soundRsamples, ~] = istft(soundRSpec, HOPSIZE, NFFT, framerate);
    soundLsamples = icwt(soundLSpec, wname, framerate);
    soundRsamples = icwt(soundRSpec, wname, framerate);
    soundSamples = [soundLsamples; soundRsamples]';
    outputFilename = sprintf('C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Output Sounds SMProto2\\%s_%d.wav', FILE, i);
    audiowrite(outputFilename, soundSamples, framerate);
end

clear i Lspectrum Rspectrum meanAmps nonzeroFreqs outputFilename soundFreqBins soundLsamples soundLsamples soundLSpec soundRsamples soundRSpec soundSamples;