function [data,fs] = read_ecg(fileName)
% Read .ecg file 
%
% Input:      
% fileName - the filename including the path
%
% Ouput:
% data     - ECG signal
% fs       - sampling frequency

% Open file for reading
fid = fopen(fileName,'r');

if ne(fid,-1) % Continue if it was possible to read the file
    
    % Magic number
    magicNumber = fread(fid, 8, 'char'); %#ok % First eight values
    
    % Get checksum
	checksum = fread(fid, 1, 'uint16'); %#ok
	
	% Read header
    Var_length_block_size = fread(fid, 1, 'long');
    ishneHeader.Sample_Size_ECG = fread(fid, 1, 'long');	
    Offset_var_lenght_block = fread(fid, 1, 'long'); %#ok
    Offset_ECG_block = fread(fid, 1, 'long');
    File_Version = fread(fid, 1, 'short'); %#ok
    First_Name = fread(fid, 40, 'char'); %#ok  		 							        								
    Last_Name = fread(fid, 40, 'char'); %#ok  									        								
    ID = fread(fid, 20, 'char'); %#ok  									        								
    Sex = fread(fid, 1, 'short'); %#ok
    Race = fread(fid, 1, 'short'); %#ok
    Birth_Date = fread(fid, 3, 'short'); %#ok	
    Record_Date = fread(fid, 3, 'short'); %#ok	
    File_Date = fread(fid, 3, 'short'); %#ok	
    Start_Time = fread(fid, 3, 'short'); %#ok	
    ishneHeader.nbLeads = fread(fid, 1, 'short');
    Lead_Spec = fread(fid, 12, 'short'); %#ok	
    Lead_Qual = fread(fid, 12, 'short'); %#ok	
    ishneHeader.Resolution = fread(fid, 12, 'short');	
    Pacemaker = fread(fid, 1, 'short'); %#ok	
    Recorder = fread(fid, 40, 'char'); %#ok
    ishneHeader.Sampling_Rate = fread(fid, 1, 'short');	
    Proprietary = fread(fid, 80, 'char'); %#ok
    Copyright = fread(fid, 80, 'char'); %#ok
    Reserved = fread(fid, 88, 'char'); %#ok
    
    % Read variable_length block
    varblock = fread(fid, Var_length_block_size, 'char'); %#ok
    
    % Get data at start
    fseek(fid, Offset_ECG_block, 'bof');
    
    % Read the signal
    numSample = ishneHeader.Sample_Size_ECG;
    data = fread(fid, [ishneHeader.nbLeads, numSample], 'int16')';
    
    % Define the sampling frequency
    fs = ishneHeader.Sampling_Rate;
     
    fclose(fid);
else
     data=[];
     fs=[];
 end
