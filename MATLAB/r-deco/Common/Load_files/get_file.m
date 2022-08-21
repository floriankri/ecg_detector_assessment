function [info, data, fs, check] = get_file(location, folder)
% Input
% location      - Location where the variable should be loaded from
% ('folder' or 'workspace')
%
% Output
% info          - Cell that contains: 
%                        - filepath   
%                        - variable (if present)
%
% data          - Selected file/variable
% fs            - If defined in the file, the sampling frequency is extracted
% check         - This is 1 when everything is executed correctly, 0 otherwise

%% Check input
if nargin == 0
    location = 'folder';
elseif nargin > 2
    errordlg('Too many input arguments!')
else
    % Set it to lower letters
    location = lower(location);
    
    % Check correctness
    if ~sum(strcmp(location,["folder","workspace"]))
        errordlg('Input arguments are not supported!')
    end
end

%% Pre-allocate
info = cell(1,2);
data = [];
fs = [];
check = 1;

%% Perform actions based on the input
try
    switch location 
        case 'folder'
            %% Get the position and name of the file
            % Define file and pathname
            
            try
                [FileName,PathName] = uigetfile({'*.MAT';'*.ECG';'*.TXT';'*.EDF';'*.CSV';'*.XML'},'Select the ECG signal',folder);
            catch
                Filename = NaN; %#ok
            end

            % If the window was not closed, check the file format
            if ~isnumeric(FileName) && ~isnumeric(PathName)

                % Get the file name
                file = fullfile(PathName,FileName);

                % Store the file path
                info{1} = file;

                % Get the format
                [~,~,format] = fileparts(FileName);

                switch lower(format)
                    case '.mat'
                        % Import data
                        data = load(file);

                        % Get the variables
                        variables = fieldnames(data);

                        if length(variables) ~= 1
                            % Select the correct variable
                            [selection,~] = listdlg('PromptString','Select the ECG file:',...
                                'ListString',variables,...
                                'SelectionMode','single');
                        else
                            selection = 1;
                        end

                        % Store the selected data variable
                        data = data.(variables{selection});

                        % Store the name of the variable
                        info{2} = variables{selection};
                        
                    case '.ecg'
                        % Import data
                        [data, fs] = read_ecg(file);
                        
                    case '.txt'
                        % Import data
                        temp = importdata(file);
                        
                        if isstruct(temp)
                            data = temp.data;
                        elseif iscell(temp)
                            % Get the nr of collumns
                            nrCol = numel(strsplit(temp{1}));            
                            data = importtextfile(file, nrCol);

                        end
                        
                    case '.csv'
                        % Import data
                        data = open(file);
                        
                        if isempty(data)
                            check = 0;
                        end
                    case '.edf'
                        % Import data
                        [data_temp, header] = read_edf(file);

                        % Sort the labels
                        [labels,idx] = sort(header.labels);

                        % Show all labels
                        [selection,~] = listdlg('PromptString','Select the ECG file(s):',...
                            'ListString',labels);

                        % Get the data
                        for ii = 1:length(selection)
                            try
                                data(:,ii) = data_temp{1,idx(selection(ii))}; %#ok
                            catch
                                data(ii,:) = data_temp{1,idx(selection(ii))}; %#ok
                            end
                        end

                        % Select the sampling frequency
                        try
                            fs = header.samplerate(idx(selection(1)));
                        catch
                        end

                        % Delete temporary data
                        clear data_temp
                    case '.h5'
                        % TODO
            %             [data,fs]=read_H5(file);
%                     case '.hea'
%                         file = listing(i);
%                         name = file.name;
% 
%                         recordName = FileName(1:end-4);
% 
%                         [ann,~,~,~,~,~] = rdann(recordName,'qrs');
%                         [ecg, fs, ~] = rdsamp(recordName);
                    case '.xml'
                        % TODO
%                         data_temp = parseXML(file);
                        
                end
            else
                check = 0;
            end

        case 'workspace'
            % Select the filepatch
            info{1} = 'Workspace';

            % Select the variables from the workspace
            variables = evalin('base','who');
            
            if ~isempty(variables)

                % Make a selection if multiple variables are present
                if length(variables) > 1
                    [selection,~] = listdlg('PromptString','Select the ECG file:',...
                        'ListString',variables,...
                        'SelectionMode','single');
                else
                    selection = 1;
                end

                % Load the variable
                data = evalin('base',variables{selection}); 

                % Store the variable name
                info{2} = variables{selection};
            else
                check = 0;
                
                % Prompt a warning dialog box
                warndlg('Workspace is empty.')
            end    
    end % switch
catch EM
    % Throw warning dlg
    warndlg(EM)
    
    check = 0;
end
