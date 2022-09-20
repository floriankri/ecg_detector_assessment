clear all

[signal,Fs,tm]=rdsamp('mitdb/100', [1 2],10000);

data(1).ekg.fs = Fs;
data(1).ekg.v = signal;

data(2).ekg.fs = Fs;
data(2).ekg.v = signal;

run_annotation(data)