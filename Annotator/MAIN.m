clear all

version = 'v005';

[signal,Fs,tm]=rdsamp('mitdb/100', 1,10000);

data(1).ekg.fs = Fs;
data(1).ekg.v = signal;

[signal,Fs,tm]=rdsamp('mitdb/102', 1,10000);

data(2).ekg.fs = Fs;
data(2).ekg.v = signal;

[signal,Fs,tm]=rdsamp('mitdb/104', 1,10000);

data(3).ekg.fs = Fs;
data(3).ekg.v = signal;


run_annotation(data,version)