clear all

version = 'v008';

currentFolder = pwd;
database_path = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\TestDatabases\annotator_test\';


cd(database_path)

name = '100';
[signal,Fs,~]=rdsamp([filesep, name,'.dat']);

data(1).ekg.fs = Fs;
data(1).ekg.v = signal;
data(1).name = name;

name = '102';
[signal,Fs,~]=rdsamp([filesep, name,'.dat']);

data(2).ekg.fs = Fs;
data(2).ekg.v = signal;
data(2).name = name;

name = '104';
[signal,Fs,~]=rdsamp([filesep, name,'.dat']);

data(3).ekg.fs = Fs;
data(3).ekg.v = signal;
data(3).name = name;

% files = dir(fullfile(database_path,['*.dat']));
% for i = 1 : length(files)
%     name = files(i).name(1:(end-4));
%     [signal,Fs,~]=rdsamp([filesep, name,'.dat']);
%     data(i).ekg.fs = Fs;
%     data(i).ekg.v = signal;
%     data(i).name = name;
% end

cd(currentFolder)
run_annotation(data,version)