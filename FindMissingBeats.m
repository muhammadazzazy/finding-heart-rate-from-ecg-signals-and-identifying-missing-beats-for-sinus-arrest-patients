function [missingBeatsTimestamps] = FindMissingBeats(file, N)
%file is the ECG of a patient with Sinus Arrest
%N is the moving average window size

% The sampling rate is 256 Hz.
fs = 256;
T = 1/fs;

% Your function could make use of the function you implemented in Part 1.
[RWaveTimestamps , RRIntervals] = DetectQRS(file, N);
average = sum(RRIntervals)/length(RRIntervals);
standardDeviation = std(RRIntervals);
count = sum(RRIntervals < average + standardDeviation);
total = sum(RRIntervals(RRIntervals < average + standardDeviation));
mean = (total/count)*10^(-3);

I = find(RRIntervals > average + standardDeviation);

numberOfMissingBeats = sum(RRIntervals > average + standardDeviation);
missingBeatsTimestamps = zeros(numberOfMissingBeats, 1);
for i = 1:1:length(missingBeatsTimestamps)
    missingBeatsTimestamps(i) = RWaveTimestamps(I(i))+mean;
end
missingBeatsTimestamps = missingBeatsTimestamps/T;
end