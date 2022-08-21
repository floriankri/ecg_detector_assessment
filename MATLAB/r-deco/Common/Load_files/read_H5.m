% Read .ecg file 
%
% Input:      
% file     - the file (name and path)
%
% Ouput:
% data     - ECG signal
% fs       - sampling frequency

function [data,fs]=read_H5(file)

% % Get the labels
% data=h5info(file,'/Stingray');
data=h5info(file);

[data,check]=loop_through_H5_structure(data);
