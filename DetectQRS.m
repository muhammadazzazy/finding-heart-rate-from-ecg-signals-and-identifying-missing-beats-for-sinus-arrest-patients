function [RWaveTimestamps, RRIntervals] = DetectQRS(file, N)
%N is the moving average window size

% The sampling rate of the ECG signal is 256 Hz
fs = 256;
T = 1/fs;
fid = fopen(file);
ECGSignal = cell2mat(textscan(fid, '%f'));
B = ECGSignal;
fclose(fid);

% A figure showing the first 1500 samples of the ECG signal before noise
% filtering.
figure(1);
n = 0:1:1499;
plot(n, ECGSignal(1:1500));
xlabel("n");
ylabel("Amplitude");
numberOfSamples = length(ECGSignal);

% A typical noisy ECG signal will appear along with a 50Hz (or 60 Hz)
% mains noise
% A notch filter with a filtered frequency of 50Hz can be used to remove
% the baseline 50Hz
% stop frequency of 50 Hz
stopFrequency = 50;
w0 = stopFrequency/(fs/2);
bw = w0/35;
[b, a] = iirnotch(w0, bw);
ECGSignal = filtfilt(b, a, ECGSignal);
% The ECG signal is further filtered using a Band-pass filter with a
% bandwidth of 0.1-45 Hz
% All frequency values are in Hz.
order = 6;  % Order
Fc1 = 0.1;  % First Cutoff Frequency
Fc2 = 45;   % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', order, Fc1, Fc2, fs);
Hd = design(h, 'butter');
ECGSignal = filter(Hd, ECGSignal);

% A figure showing the first 1500 samples of the ECG signal after noise
% filtering.
figure(2);
plot(n, ECGSignal(1:1500));
xlabel("n");
ylabel("Amplitude");

% One way to overcome baseline drifts is to differentiate the ECG signal
% since we are looking for large slopes in the signal
x = ECGSignal;
y = zeros(numberOfSamples, 1);
y(1) = (1/(8*T))*(2*x(2)+x(3));
y(2) = (1/(8*T))*(-2*x(1)+2*x(3)+x(4));
for n = 3:1:numberOfSamples-2
    y(n) = (1/(8*T))*(-x(n-2)-2*x(n-1)+2*x(n+1)+x(n+2));
end
y(numberOfSamples-1) = (1/(8*T))*(-2*x(numberOfSamples-2)-x(numberOfSamples-3)+2*x(numberOfSamples));
y(numberOfSamples) = (1/(8*T))*(-2*x(numberOfSamples-1)-x(numberOfSamples-2));

% The next step is to square the derivative. This makes all data points
% positive and does nonlinear amplification of the output of the derivative
% emphasizing the higher frequencies (i.e., predominantly the ECG
% frequencies)
x = y;
for n = 1:1:numberOfSamples
    y(n) = (x(n))^2;
end

% The next step is to smooth the squared signal using a moving-average
% window
% y(nT) = (1/N)[x(nT-(N-1)T)+x(nT-(N-2)T)
% +...+x(nT)]
% where N is the number of samples in the width of the moving average
% window
% N should be approximately the same as the widest possible QRS complex
x = y;
y = zeros(numberOfSamples, 1);
for n = 1:1:numberOfSamples
    for i = 1:1:N
        if(n >= N)
            y(n) = y(n) + x(n-N+i);
        else
            y(n) = x(n);
        end
    end
    y(n) = y(n)/N;
end
ECGSignal = y;

% Zero padding
for i = 1:1:fs/2
    ECGSignal(i) = 0;
end
for i = numberOfSamples-fs/2:1:numberOfSamples
    ECGSignal(i) = 0;
end

% The final step is to set a threshold on the moving average output. The
% peak value above the threshold within the moving average window is
% approximately the R wave
A = zeros(fs/2, 1);
for i = fs:1:fs*2
    if(floor(ECGSignal(i)) <= 190)
        A(i) = ECGSignal(i);
    end
end
threshold = 15*std(A);
[pks, locs] = findpeaks(ECGSignal, "MinPeakHeight", threshold, "MinPeakProminence", N);
RWaveTimestamps = zeros(length(pks), 1);
for i = 1:1:length(RWaveTimestamps)
    RWaveTimestamps(i) = locs(i)*T;
end

% The sequence of RR intervals - that is, all intervals between adjacent QRS
% complexes resulting from sinus node depolarizations - forms the RR
% interval time series or RR tachogram
RRIntervals = zeros(length(pks)-1, 1);
for i = 2:1:length(RWaveTimestamps)
    RRIntervals(i-1) = (RWaveTimestamps(i)-RWaveTimestamps(i-1))*10^3;
end

% A corresponding sequence of instantaneous heart rate is defined as ff_i =
% 1/RR_i
heartRate = zeros(length(RRIntervals), 1);
for i = 1:1:length(heartRate)
    heartRate(i) = 1/RRIntervals(i);
end

% A figure showing the first 1500 samples of the ECG signal with an "*"
% marking the detected R waves.
n = 0:1:1499;
figure(3);
plot(n, ECGSignal(1:1500));
hold on
plot(n(locs(1:5)), ECGSignal(locs(1:5)), "*");
xlabel("n");
ylabel("Amplitude");

% A figure showing the first 1500 samples of the ECG signal with an "*" marking the detected
% R waves for N = 25 but without noise filtering.
figure(4);
plot(n, B(1:1500));
hold on
plot(n(locs(1:5)), B(locs(1:5)), "*");

% A plot of the RR intervals with Beat number on the x-axis and RR interval
% in msec on the y-axis in the case of N = 25.
figure(5);
x = 0:1:length(RRIntervals)-1;
plot(x, RRIntervals);
xlabel("Beat number");
ylabel("RR interval (msec)");
end