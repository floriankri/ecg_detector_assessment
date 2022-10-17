clear all

version = 'v007';

path = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\TestDatabases\annotator_test\';

cd(path)
name = '100';
[signal,Fs,tm]=rdsamp([filesep, name,'.dat'], 1,10000);

data(1).ekg.fs = Fs;
data(1).ekg.v = signal;
data(1).name = name;

name = '102';
[signal,Fs,tm]=rdsamp([filesep, name,'.dat'], 1,10000);

data(2).ekg.fs = Fs;
data(2).ekg.v = signal;
data(2).name = name;

name = '104';
[signal,Fs,tm]=rdsamp([filesep, name,'.dat'], 1,10000);

data(3).ekg.fs = Fs;
data(3).ekg.v = signal;
data(3).name = name;

cd('C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\Annotator')
run_annotation(data,version)


% [filename,directoryname] = uigetfile('*.hea;*.edf;*.mat','Select signal header file:');
% cd(directoryname)
% fext=filename(end-3:end);
% tmp=dir(['*' fext]);
% N=length(tmp);
% records=cell(N,1);
% for n=1:N
%     fname=tmp(n).name;
%     if(strcmp(fext,'.edf'))
%         records(n)={fname};
%     else
%         records(n)={fname(1:end-4)};
%     end
%     if(strcmp(fname,filename))
%         current_tmp=n;
%     end
% end