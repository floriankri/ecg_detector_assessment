function varargout = R_peak_detection(varargin)
% R_PEAK_DETECTION MATLAB code for R_peak_detection.fig
%      R_PEAK_DETECTION, by itself, creates a new R_PEAK_DETECTION or raises the existing
%      singleton*.
%
%      H = R_PEAK_DETECTION returns the handle to a new R_PEAK_DETECTION or the handle to
%      the existing singleton*.
%
%      R_PEAK_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in R_PEAK_DETECTION.M with the given input arguments.
%
%      R_PEAK_DETECTION('Property','Value',...) creates a new R_PEAK_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before R_peak_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to R_peak_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 06/05/2019   JM      Initial version
%
% Copyright (c) 2019,  Jonathan Moeyersons, KULeuven-ESAT-STADIUS
%
% This software is made available for non commercial research purposes only
% under the GNU General Public License. However, notwithstanding any
% provision of the GNU General Public License, this software may not be
% used for commercial purposes without explicit written permission after
% contacting jonathan.moeyersons@esat.kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Edit the above text to modify the response to help R_peak_detection

% Last Modified by GUIDE v2.5 14-May-2019 11:53:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @R_peak_detection_OpeningFcn, ...
    'gui_OutputFcn',  @R_peak_detection_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before R_peak_detection is made visible.
function R_peak_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to R_peak_detection (see VARARGIN)

% Choose default command line output for R_peak_detection
handles.output = hObject;

% Set figure
imshow('R_peak_example.png','Parent',handles.axes_Rpeak_example)

% Empty the folder variable
handles.data.folder = [];

% Set default variables
handles = reset(hObject, eventdata, handles);

% Get the data if it is given
if ~isempty(varargin) && (mod(length(varargin),6) == 0 || mod(length(varargin),4) == 0)
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'signal'
                data = varargin{ii+1};
            case 'fs'
                fs = varargin{ii+1};
            case 'Rpeaks'
                Rpeaks = varargin{ii+1};
        end
    end
    
    % Put the data in collumns
    [R_temp , K_temp] = size(data);
    if R_temp < K_temp
        data = data';
        Length = K_temp;
        channels = R_temp;
    else
        Length = R_temp;
        channels = K_temp;
    end
    
    % Compute the time
    time = seconds((1:Length)/fs);
    
    % Get the length of the signal
    temp_end = time(end);
    temp_end.Format = 'hh:mm:ss';
    
    % Set the data parameters
    handles.data.folder = cd;
    handles.data.filepath = 'Workspace';
    handles.data.variable = 'Signal';
    handles.data.signal.original = data;
    handles.data.signal.filtered = data;
    handles.data.fs = fs;
    handles.data.duration_recording = temp_end;
    handles.data.start_analysis = duration([0 0 0]);
    handles.data.duration_analysis = temp_end;
    handles.data.range = temp_end;
    handles.data.channels = channels;
    
    % Get the R-peak locations
    if mod(length(varargin),6) == 0
        handles.data.R{1} = Rpeaks;
    end
    
    % Update the gui
    handles = update_gui(handles);
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
    
    % Add the listener of the slider
    handles.listener = addlistener(handles.slider_range ,'Value','PostSet',...
        @Slide);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes R_peak_detection wait for user response (see UIRESUME)
% uiwait(handles.figure_main);


% --- Outputs from this function are returned to the command line.
function varargout = R_peak_detection_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_reset.
function handles = pushbutton_reset_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if a signal has been loaded and the results are not saved yet
choice = 'Yes';
if ~handles.check.save && ~isempty(handles.data.signal.original)
    choice = questdlg('Current analysis has not been saved. Are you sure you want to reset?', ...
        'Reset', ...
        'Yes','No','Cancel','Yes');
end

% Next step depends on the answer of the previous question
switch choice
    case 'Yes'
        % Reset the gui
        handles = reset(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);


function handles = reset(hObject, eventdata, handles)%#ok
% hObject    handle to selected uicontrol
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset the menu
set(handles.menu_save_session,'enable','off')
set(handles.menu_export_results,'enable','off')
set(handles.menu_save_results,'enable','off')

% Reset the zoom
zoom off

% Reset the windowbuttons
set(handles.figure_main,'WindowButtonMotionFcn','',...
    'WindowButtonDownFcn','',...
    'WindowButtonUpFcn','');

% Reset the resize flag
handles.uiLimitSet = false;

% Reset the data panel
set(handles.text_filepath,'string','')
set(handles.text_variable,'string','')
set(handles.text_sampling_frequency,'string','')
set(handles.text_signal_duration,'string','')
set(handles.text_channels,'string','')
set(findall(handles.uipanel_data,'-property','enable'),'enable','on')
set(handles.pushbutton_reset,'enable','off')

% Reset the filter panel
set(findall(handles.uipanel_filter,'Style','edit'),'enable','off',...
    'string','')
set(findall(handles.uipanel_filter,'Style','pushbutton'),'enable','off')
set(handles.checkbox_show_original,'Value',0,...
    'enable','off')

% Reset the analysis period panel
set(handles.edit_start,'string','00:00:00');
set(handles.edit_duration,'string','00:00:00');

set(handles.togglebutton_analysis_period,'enable','off')
set(handles.edit_start,'enable','off')
set(handles.edit_duration,'enable','off')

% Reset the R-peak detection panel
set(findall(handles.uipanel_R_peak_detection,'-property','enable'),'enable','off')

% Reset the parameters
handles.data.parameters.Welch{1} = []; % Window width
handles.data.parameters.Welch{2} = []; % Overlap
handles.data.parameters.Welch{3} = []; % Nfft

handles.data.parameters.Filter{1} = 2; % High pass order
handles.data.parameters.Filter{2} = 4; % Low pass order

handles.data.parameters.Rpeak{1} = 300;
handles.data.parameters.Rpeak{2} = 100;
handles.data.parameters.Rpeak{3} = 1;
handles.data.parameters.Rpeak{4} = 0;
handles.data.parameters.Rpeak{5} = 0;

% Reset the R-peak correction panel
set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)
set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','off')

% Reset the range
set(handles.edit_range,'string','00:01:00',...
    'enable','off')
set(handles.pushbutton_plus,'enable','off')
set(handles.pushbutton_minus,'enable','off')
set(handles.text_range,'enable','off')

% Reset the slider
try
    % Delete listener
    delete(handles.listener)
catch
end
set(handles.slider_range,'value',0,...
    'string','normal',...
    'enable','off')

% Reset the channels listbox
set(handles.listbox_channels,'string','All',...
    'value',1,...
    'enable','off')

% Reset the tachogram units
set(handles.popupmenu_units_tachogram,'value',1,...
    'enable','off')

% Disable the fixed y-limits check box
set(handles.checkbox_fix_ylimits,'value',0,...
    'enable','off')

% Clear the data and results variables
[handles.data.filepath,handles.data.variable,handles.data.signal.original,...
    handles.data.signal.filtered,handles.data.R,handles.data.fs,...
    handles.data.duration_recording,...
    handles.data.channels,handles.data.hp,handles.data.lp,...
    handles.data.stop] = deal([]);

% Reset the analysis window
handles.data.start_analysis = seconds(0);
handles.data.duration_analysis = seconds(60);
handles.data.range = seconds(60);

handles.data.units = 1;

[handles.graph.signal.original,handles.graph.signal.filtered,...
    handles.graph.R,handles.graph.RR] = deal([]);

% (Re)set the check
handles.check.save = 0;

% (Re)set the channel labels
handles.labels = {'Channel_1','Channel_2','Channel_3','Channel_4',...
    'Channel_5','Channel_6','Channel_7','Channel_8','Channel_9','Channel_10',...
    'Channel_11','Channel_12','Channel_13'};

% Reset to default colors
handles.color = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840;...
    rand(6,3)];

% Clear all axes
cla(handles.axes_ECG)
cla(handles.axes_tachogram)

% Plot random something on the axes
plot(handles.axes_ECG,seconds([0 60]),[0 1],...
    'visible','off');

plot(handles.axes_tachogram,seconds([0 60]),[0 1],...
    'visible','off');

% Adjust the xtickformat
xtickformat([handles.axes_ECG handles.axes_tachogram],'hh:mm:ss')

% Adjust the x- and y-limits
xlim([handles.axes_ECG handles.axes_tachogram],seconds([0 60]))
ylim([handles.axes_ECG handles.axes_tachogram],[0 1])

% Linkaxes
linkaxes([handles.axes_ECG,handles.axes_tachogram],'x')

% Define the axis labels
xlabel(handles.axes_tachogram,{'Time (hh:mm:ss)'})
ylabel(handles.axes_ECG,'ECG (a.u.)')
ylabel(handles.axes_tachogram,'RR (ms)')

% Reset the context menu and remove the gridlines
set(handles.context_x_grid,'Checked','off')
set(handles.context_y_grid,'Checked','off')
set(handles.context_add,'Enable','off')
set(handles.axes_ECG,'XGrid','off')
set(handles.axes_ECG,'YGrid','off')
set(handles.axes_tachogram,'XGrid','off')
set(handles.axes_tachogram,'YGrid','off')

% Clear all axes again
cla(handles.axes_ECG)
cla(handles.axes_tachogram)

% Enable the picture pushbutton
set(handles.pushbutton_picture,'enable','on')

% Set zoom and zoom postcallback
handles.uitoggletool_zoom_in.Enable = 'on';
handles.uitoggletool_zoom_out.Enable = 'on';

handles.zoom = zoom;
handles.zoom.ActionPostCallback = @zoom_postcallback;
handles.zoom.Enable = 'off';

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_load_session_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_load_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if a signal has been loaded and if the results are already saved
choice = 'Yes';
if ~handles.check.save && ~isempty(handles.data.signal.original)
    choice = questdlg('Current analysis has not been saved. Are you sure you want to load a session?', ...
        'Load session', ...
        'Yes','No','Cancel','Yes');
end

% Next step depends on the answer of the previous question
switch choice
    case 'Yes'
        % Change the pointer
        set(handles.figure_main,'Pointer','watch')
        drawnow;
        
        % Define file and pathname
        [FileName,PathName] = uigetfile('*.mat','Select the session');
        
        % If the window was not closed, check the file format
        if ~isnumeric(FileName) && ~isnumeric(PathName)
            
            % Reset everything
            handles = reset(hObject, eventdata, handles);
            
            % Get the complete file
            file = fullfile(PathName,FileName);
            
            % Import data
            temp = load('-mat',file);
            
            % Set the data parameters
            handles.data = temp.data;
            
            % Update the gui
            handles = update_gui(handles);
            
            % Adjust the y-limits
            handles = adjust_y_limits(handles);
            
            % Add the listener of the slider
            handles.listener = addlistener(handles.slider_range ,'Value','PostSet',...
                @Slide);
        end
        % Change the pointer
        set(handles.figure_main,'Pointer','arrow')
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_save_session_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_save_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Define file and pathname
[FileName,PathName] = uiputfile('.mat','Save as');

% If the window was not closed, save the results
if ~isnumeric(FileName) && ~isnumeric(PathName)
    data = handles.data;
    
    % Adjust R-peaks if one of the adjustment buttons is still selected
    for ii = 1:handles.data.channels
        try
            data.R{ii} = handles.graph.R.(handles.labels{ii}).XData;
        catch
            data.R{ii} = [];
        end
    end
    
    % Save session
    N = strcat(PathName,FileName);
    save(N,'data')
end


% --------------------------------------------------------------------
function menu_save_results_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_save_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Define file and pathname
[FileName,PathName] = uiputfile('.mat','Save as');

% If the window was not closed, save the results
if ~isnumeric(FileName) && ~isnumeric(PathName)
    
    % Loop over the channels
    for ii = 1:handles.data.channels
        
        % Store the intervals and locations in cells
        try
            temp = handles.graph.R.(handles.labels{ii}).XData;
            temp.Format = 'hh:mm:ss.SSSS';
            data.R_loc{ii}  = temp;
            data.RR_int{ii} = handles.graph.RR.(handles.labels{ii}).YData;
        catch
            data.R_loc{ii}  = [];
            data.RR_int{ii} = [];
        end
    end
    
    % Save results
    N = strcat(PathName,FileName);
    save(N,'data')
end


% --------------------------------------------------------------------
function menu_to_excel_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_to_excel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Creat an excel sheet as output
create_excel_output(handles)

% Sampling frequency
% Total duration recording
% # Channels
% Total duration analysis period
% Total nr of beats of recording
% Total nr of beats of analysis
% Mean HR
% Median HR
% Min HR
% Max HR


% --------------------------------------------------------------------
function menu_to_workspace_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_to_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Arrange output
for ii = 1:handles.data.channels
    data.(handles.labels{ii}).signal = handles.graph.signal.filtered.(handles.labels{ii}).YData;
    data.(handles.labels{ii}).R = round(handles.data.fs*seconds(handles.graph.R.(handles.labels{ii}).XData-handles.data.start_analysis));
    data.(handles.labels{ii}).RR = handles.graph.RR.(handles.labels{ii}).YData;
end

% Send to workspace
assignin('base','data',data);

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_quit_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if a signal has been loaded and if the results are already saved
choice = 'Yes';
if ~handles.check.save && ~isempty(handles.data.signal.original)
    choice = questdlg('Current analysis has not been saved. Are you sure you want to quit?', ...
        'Quit', ...
        'Yes','No','Cancel','Yes');
end

% Next step depends on the answer of the previous question
switch choice
    case 'Yes'
        delete(handles.figure_main)
end


% --------------------------------------------------------------------
function menu_preferences_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to menu_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~verLessThan('matlab','9.4')
    % Get all to be disabled objects
    hObjects = findall(handles.figure_main,'-property','enable');
    
    % Get the state of all objects
    Objects_state = get(hObjects,'enable');
    
    % Disable everything
    set(hObjects,'enable','off');
    
    % Run the preferences app
    hPref = Preferences(handles.data.parameters);
    
    % Get the figure
    handles.preferences.Fig = hPref.PreferencesUIFigure;
    
    % Update handles structure (necesarry when figure is deleted)
    guidata(hObject, handles);
    
    % Wait for the preferences to be adjusted
    uiwait(handles.preferences.Fig)
    
    % Extract the parameters
    if isvalid(handles.preferences.Fig)
        handles.data.parameters = hPref.Parameters;
        
        % Delete the app
        close(handles.preferences.Fig)
    end
    
    % Update GUI (if the GUI still exists)
    if isvalid(handles.figure_main)
        for ii = 1:length(hObjects)
            set(hObjects(ii),'enable',Objects_state{ii})
        end
        
        % Update handles structure
        guidata(hObject, handles);
    end
else
    % Run the preferences GUI
    temp = Preferences_GUIDE(handles.data.parameters);
    
    if ~isempty(temp)
        % Store the parameters
        handles.data.parameters = temp;
        
        % Update handles structure
        guidata(hObject, handles);
    end
end


% --- Executes on button press in pushbutton_from_file.
function pushbutton_from_file_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Check if a signal has been loaded and if the results are already saved
choice = 'Yes';
if ~handles.check.save && ~isempty(handles.data.signal.original)
    choice = questdlg('Current analysis has not been saved. Are you sure you want to start a new session?', ...
        'Start new session', ...
        'Yes','No','Cancel','Yes');
end

% Next step depends on the answer of the previous question
switch choice
    case 'Yes'
        % Reset the manual corrections
        set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)
        
        % Execute the radiobutton callback
        radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
        
        % Select file
        if strcmp(hObject.Tag(end-3:end),'file')
            [info, data, fs, check] = get_file('folder',handles.data.folder);
        else
            [info, data, fs, check] = get_file('workspace');
        end
        
        if check
            % If the selected file is a structure, dig deeper
            if isa(data,'struct')
                if length(fieldnames(data)) > 1 || isa(data.(char(fieldnames(data))),'struct')
                    [data , check] = loop_through_structure(data);
                else
                    data = data.(char(fieldnames(data)));
                end
            end
            
            if check
                % Check if the selected file is a tensor
                if ~(length(size(data)) > 2)
                    
                    % Check if the selected file is large enough
                    if max(size(data)) > 500
                        
                        % Check if the file contains nans
                        if sum(isnan(data(:))) == 0
                            
                            % If the sampling frequency has not been defined
                            if isempty(fs)
                                [fs , check] = get_sampling_frequency;
                            end
                            
                            if check
                                % Set the signal in the correct form (collumn vectors)
                                [R_temp , K_temp] = size(data);
                                if R_temp < K_temp
                                    data = data';
                                    Length = K_temp;
                                    channels = R_temp;
                                else
                                    Length = R_temp;
                                    channels = K_temp;
                                end
                                
                                % Compute the time
                                time = seconds((1:Length)/fs);
                                
                                % Get the length of the signal
                                temp_end = time(end);
                                temp_end.Format = 'hh:mm:ss';
                                
                                % Reset everything
                                handles = reset(hObject, eventdata, handles);
                                
                                % Set the data parameters
                                [handles.data.folder, ~, ~] = fileparts(info{1});
                                handles.data.filepath = info{1};
                                handles.data.variable = info{2};
                                handles.data.signal.original = data;
                                handles.data.signal.filtered = data;
                                handles.data.fs = fs;
                                handles.data.duration_recording = temp_end;
                                handles.data.start_analysis = duration([0 0 0]);
                                handles.data.duration_analysis = temp_end;
                                handles.data.range = temp_end;
                                handles.data.channels = channels;
                                
                                % Update the gui
                                handles = update_gui(handles);
                                
                                % Adjust the y-limits
                                handles = adjust_y_limits(handles);
                                
                                % Add the listener of the slider
                                handles.listener = addlistener(handles.slider_range ,'Value','PostSet',...
                                    @Slide);
                            end
                        else
                            errordlg('The selected input file contains NaNs','File Error','modal')
                        end
                    else
                        errordlg('The selected input file is too small','File Error','modal')
                    end
                else
                    errordlg('The selected input file is not supported','File Error','modal')
                end
            end
        end
end

% Change the pointer
set(handles.figure_main,'Pointer','arrow')

% Update handles structure
guidata(hObject, handles);


function edit_high_pass_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_high_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_high_pass as text
%        str2double(get(hObject,'String')) returns contents of edit_high_pass as a double


% --- Executes during object creation, after setting all properties.
function edit_high_pass_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_high_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_low_pass_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_low_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_low_pass as text
%        str2double(get(hObject,'String')) returns contents of edit_low_pass as a double


% --- Executes during object creation, after setting all properties.
function edit_low_pass_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_low_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_stop_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stop as text
%        str2double(get(hObject,'String')) returns contents of edit_stop as a double


% --- Executes during object creation, after setting all properties.
function edit_stop_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_check_psd.
function pushbutton_check_psd_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_check_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Define the looping variable
if listbox_value == 0
    loop_var = 1:handles.data.channels;
else % Only the selection
    loop_var = listbox_value;
end

% Get the start and stop of the analysis window in samples
start = max([1 round(handles.data.fs*seconds(handles.data.start_analysis))]);
stop = start + min([length(handles.data.signal.filtered) round(handles.data.fs*seconds(handles.data.duration_analysis))])-1;

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Create figure
hfig = figure('Name','Power Spectrum','NumberTitle','off');

% Create tabgroup
htgroup = uitabgroup('Parent', hfig,'TabLocation', 'top');

% Remove the menu bar
hfig.MenuBar = 'none';

% Adjust the toolbar
hfig.ToolBar = 'figure';

% Get the toolbar
htb = findall(hfig,'Type','uitoolbar');

% Find all tools
tools = findall(htb);

% Check for the version
if verLessThan('matlab', '9.5')
    % Delete unnecessary tools
    delete(tools([2:7 9:10 13:17]))
    
    % Remove the separator
    set(tools(12),'Separator','off')
else
    % Delete all tools
    delete(tools)
end

% Loop over the selected channels
for ii = loop_var
    
    % Create panel
    h.(strcat('panel',num2str(ii))) = uipanel('Position',[0 0 1 1]);
    
    % Create tab
    h.(strcat('tab',num2str(ii))) = uitab('Parent', htgroup, 'Title', strcat("Channel ",num2str(ii)));
    
    % Set parent
    set(h.(strcat('panel',num2str(ii))),'Parent',h.(strcat('tab',num2str(ii))))
    
    % Create axes
    h.(strcat('axes',num2str(ii))) = axes('Parent',h.(strcat('panel',num2str(ii))));
    
    % Define signal and original signal
    signal = handles.data.signal.filtered(start:stop,ii);
    original = handles.data.signal.original(start:stop,ii);
    
    % Compute the power spectral density with the p-welch method
    window_width = round(handles.data.parameters.Welch{1}*handles.data.fs);
    
    if window_width >= length(signal)
        % Give a warning and use the default values
        warndlg({'Window width is larger than the signal.','Default values are used.'})
        
        % Adjust parameters
        window_width = [];
        noverlap = [];
        nfft = [];
    else
        % Get other parameters
        noverlap = round(window_width/100*handles.data.parameters.Welch{2});
        nfft = round(handles.data.parameters.Welch{3}*handles.data.fs);
    end
    
    % Compute the power spectral density with the p-welch method
    [pxx,f] = pwelch(signal,window_width,noverlap,nfft,handles.data.fs);
    
    % Make the plot
    plot(f,10*log10(pxx),'linewidth',1.5,...
        'color',handles.color(ii,:))
    
    % Compute the power spectral density with the p-welch method
    [pxx,f] = pwelch(original,window_width,noverlap,nfft,handles.data.fs);
    
    % Make the plot
    hold on;
    plot(f,10*log10(pxx),'linewidth',1.5,...
        'color',[handles.color(ii,:),0.2])
    
    % Make a legend
    legend('Filtered','Original','Location','northeast')
    
    % Set the x-limits
    xlim([0 f(end)])
    
    % Set the labels
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
end

% % Change windowstyle
% set(hfig,'WindowStyle','modal')

% Change the pointer
set(handles.figure_main,'Pointer','arrow')
drawnow;


% --- Executes on button press in pushbutton_filter.
function pushbutton_filter_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if the R-peaks have been detected already
choice = 'Yes';
if ~isempty(handles.data.R)
    choice = questdlg('This action will remove the detected R-peaks. Are you sure you want to continue?', ...
        'Filter signal', ...
        'Yes','No','Cancel','Yes');
end

switch choice
    case 'Yes'
        
        % Get all strings
        hp_string = get(handles.edit_high_pass,'string');
        lp_string = get(handles.edit_low_pass,'string');
        stop_string = get(handles.edit_stop,'string');
        
        % Check if anything is present
        if ~isempty(hp_string) || ~isempty(lp_string) || ~isempty(stop_string)
            % Replace comma's
            hp_string = replace_comma(hp_string);
            lp_string = replace_comma(lp_string);
            stop_string = replace_comma(stop_string);
            
            % Get all values
            hp = str2double(hp_string);
            lp = str2double(lp_string);
            stop = str2double(stop_string);
            
            % Check if it is possible
            if (isempty(hp_string) || (~isnan(hp) && hp < handles.data.fs/2 && hp > 0)) &&...
                    (isempty(lp_string) || (~isnan(lp) && lp < handles.data.fs/2 && lp > 0)) &&...
                    (isempty(stop_string) || (~isnan(stop) && stop < handles.data.fs/2 && stop > 0))
                
                % Change the pointer
                set(handles.figure_main,'Pointer','watch')
                drawnow;
                
                % Reset the manual corrections
                set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)
                
                % Execute the radiobutton callback
                radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
                
                % Filter the signal
                [signal,check] = filtersignal(handles.data.signal.original,handles.data.fs,[~isnan(hp),~isnan(lp),~isnan(stop),0],handles.data.parameters.Filter,hp,lp,stop);
                
                if check
                    % Adjust the model
                    handles.data.signal.filtered = signal;
                    
                    if ~isnan(hp)
                        handles.data.hp     = hp;
                        set(handles.edit_high_pass,'string',hp)
                    end
                    if ~isnan(lp)
                        handles.data.lp     = lp;
                        set(handles.edit_low_pass,'string',lp)
                    end
                    if ~isnan(stop)
                        handles.data.stop     = stop;
                        set(handles.edit_stop,'string',stop)
                    end
                    
                    % Remove previously detected R-peaks if present
                    if ~isempty(handles.data.R)
                        % Clear the model
                        handles.data.R = [];
                    end
                    
                    % Adjust the plots
                    handles = update_plots(handles);
                    
                    % Adjust the y-limits
                    handles = adjust_y_limits(handles);
                    
                    % Set the filter panel
                    set(handles.pushbutton_remove_filter,'enable','on')
                    set(handles.checkbox_show_original,'enable','on')
                    
                    % Set the R-peak correction panel
                    set(findall(handles.uipanel_R_peak_correction, '-property', 'enable'),'enable','off')
                    
                    % Disable the addition of R-peaks via the context menu
                    set(handles.context_add,'Enable','off')
                    
                    % Reset the menu
                    set(handles.menu_export_results,'enable','off')
                    set(handles.menu_save_results,'enable','off')
                end
                % Change the pointer
                set(handles.figure_main,'Pointer','arrow')
                
            else
                % Change the string
                set(handles.edit_high_pass,'string','')
                set(handles.edit_low_pass,'string','')
                set(handles.edit_stop,'string','')
                
                % Throw an error dialog
                errordlg('The input format is not supported','Format Error','modal')
            end
        else
            % Prompt a warning dialog box
            warndlg('Thresholds are empty.')
        end
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in pushbutton_remove_filter.
function pushbutton_remove_filter_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_remove_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if the R-peaks have been detected already
choice = 'Yes';
if ~isempty(handles.data.R)
    choice = questdlg('This action will remove the detected R-peaks. Are you sure you want to continue?', ...
        'Remove filter', ...
        'Yes','No','Cancel','Yes');
end

switch choice
    case 'Yes'
        % Reset the manual corrections
        set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)
        
        % Execute the radiobutton callback
        radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
        
        % Reset signal
        handles.data.signal.filtered = handles.data.signal.original;
        
        % Clear the R-peaks
        if ~isempty(handles.data.R)
            handles.data.R = [];
        end
        
        % Adjust the original checkbox
        set(handles.checkbox_show_original,'enable','off',...
            'Value',0)
        
        % Adjust the plots
        handles = update_plots(handles);
        
        % Adjust the y-limits
        handles = adjust_y_limits(handles);
        
        % Set the filter panel
        set(handles.pushbutton_remove_filter,'enable','off')
        
        % Set the R-peak correction panel
        set(findall(handles.uipanel_R_peak_correction, '-property', 'enable'),'enable','off')
        
        % Disable the addition of R-peaks via the context menu
        set(handles.context_add,'Enable','off')
        
        % Reset the menu
        set(handles.menu_export_results,'enable','off')
        set(handles.menu_save_results,'enable','off')
        
        % Change the string
        set(handles.edit_high_pass,'string','')
        set(handles.edit_low_pass,'string','')
        set(handles.edit_stop,'string','')
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_show_original.
function checkbox_show_original_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to checkbox_show_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_show_original

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Define the looping variable
if listbox_value == 0
    loop_var = 1:handles.data.channels;
else % Only the selection
    loop_var = listbox_value;
end

% Make the selected channels (in)visible
for ii = loop_var
    set(handles.graph.signal.original.(handles.labels{ii}),'visible',get(hObject,'Value'))
end

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Update structure
guidata(hObject,handles);


% --- Executes on button press in togglebutton_analysis_period.
function togglebutton_analysis_period_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to togglebutton_analysis_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_analysis_period


% Get the togglebutton value
togglebutton_value = get(handles.togglebutton_analysis_period,'value');

switch togglebutton_value
    case 1
        % Check if the R-peaks have been detected already
        choice = 'Yes';
        if ~isempty(handles.data.R)
            choice = questdlg('Do you want to remove the currently detected R-peaks?', ...
                'Define analysis window', ...
                'Yes','No','Cancel','Yes');
        end
        
        if ~strcmp(choice,'')
            
            % Reset the manual corrections
            set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)
            
            % Execute the radiobutton callback
            handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
            
            % Reset start analysis period
            set(handles.edit_start,'string','00:00:00')
            handles.data.start_analysis = duration([0 0 0]);
            
            % Reset duration analysis period
            set(handles.edit_duration,'string',char(handles.data.duration_recording))
            handles.data.duration_analysis = handles.data.duration_recording;
            
            % Delete the R-peaks if wanted
            if strcmp(choice,'Yes')
                handles.data.R = [];
                
                % Disable the addition of R-peaks via the context menu
                set(handles.context_add,'Enable','off')
            end
            
            % Adjust the plots
            handles = update_plots(handles);
            
            % Update structure
            guidata(hObject,handles);
            
            % Change the togglebutton text
            set(handles.togglebutton_analysis_period,'String','Apply changes')
            
            % Set slider
            if ~strcmp(char(handles.data.duration_analysis),char(handles.data.range))
                % Get the x-limits
                Xlim = xlim(handles.axes_ECG);
                
                % Get a new maximum for the slider
                new_max = round(handles.data.fs*seconds(handles.data.duration_analysis-handles.data.range));
                
                % Get the new step size
                n = round(handles.data.fs*seconds(handles.data.range))/new_max;
                
                % Adjust the slider
                set(handles.slider_range,'max',new_max,...
                    'sliderstep',[min([1 n]) min([1 2*n])],...
                    'value',round(handles.data.fs*seconds(Xlim(1)-handles.data.start_analysis)),...
                    'string','zoom',...
                    'enable','on')
                
                set(handles.slider_range,'string','normal')
            end
            
            % Create the analysis period patch
            Ylim = ylim(handles.axes_ECG);
            
            handles.temp.patch = patch(handles.axes_ECG,[0 0 0 0],...
                [Ylim(1) Ylim(1) Ylim(2) Ylim(2)],...
                'b',...
                'facealpha',0.1,...
                'visible','off');
            
            % Create the starting line
            handles.temp.start = line(handles.axes_ECG,seconds([0 0]),Ylim,...
                'visible','off');
            
            % Create the stop line
            handles.temp.stop = line(handles.axes_ECG,[handles.data.duration_recording handles.data.duration_recording],Ylim,...
                'visible','off');
            
            % Disable the data and R-peaks uipanel
            set(findall(handles.uipanel_data, '-property', 'enable'),'enable','off')
            set(findall(handles.uipanel_filter, '-property', 'enable'),'enable','off')
            set(findall(handles.uipanel_R_peak_detection, '-property', 'enable'),'enable','off')
            set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','off')
            
            % Set the window buttonmotion, down and up function
            set(handles.figure_main,'WindowButtonMotionFcn',{@drag_start_line,handles},...
                'WindowButtonDownFcn',{@define_start_line,handles},...
                'WindowButtonUpFcn',{@define_stop_line,handles});
        else
            % Reset the togglebutton
            set(handles.togglebutton_analysis_period,'value',0)
        end
    case 0
        % Get the start time
        if strcmp(handles.temp.start.Visible,'on')
            start = handles.temp.start.XData(1);
        else
            start = handles.data.start_analysis;
        end
        
        % Get the stop time
        if strcmp(handles.temp.stop.Visible,'on')
            stop = handles.temp.stop.XData(1);
        else
            stop = handles.data.duration_analysis;
        end
        
        % Adjust the analysis start
        handles.data.start_analysis = start;
        
        % Adjust the analysis duration
        handles.data.duration_analysis = stop-start;
        
        % Adjust the plots
        handles = update_plots(handles);
        
        % Adjust the x-limits
        xlim(handles.axes_ECG,[start stop])
        
        %         % Get the fix value
        %         fix = get(handles.checkbox_fix_ylimits,'Value');
        
        % Adjust the y-limits
        handles = adjust_y_limits(handles);
        
        % Adjust the range string
        set(handles.edit_range,'string',char(stop-start))
        handles.data.range = stop-start;
        
        % Disable the slider
        set(handles.slider_range,'enable','off')
        
        % Change the togglebutton text
        set(handles.togglebutton_analysis_period,'String','Define analysis period')
        
        % Disable the analysis period edit boxes
        set(findall(handles.uipanel_analysis_period,'style','edit'),'enable','off')
        
        % Enable the data and R-peaks pushbuttons, the listbox and the toggle
        % tools
        set(findall(handles.uipanel_data, '-property', 'enable'),'enable','on')
        set(findall(handles.uipanel_filter, '-property', 'enable'),'enable','on')
        
        if sum((handles.data.signal.original(:) == handles.data.signal.filtered(:))-1) == 0
            % Set the filter panel
            set(handles.pushbutton_remove_filter,'enable','off')
            set(handles.checkbox_show_original,'enable','off')
        end
        
        set(handles.pushbutton_detect_peaks,'enable','on')
        set(handles.pushbutton_load_peaks,'enable','on')
        
        if ~isempty(handles.data.R)
            set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','on')
        end
        
        zoom off
        
        % Set the window buttonmotion, down and up function
        set(handles.figure_main,'WindowButtonMotionFcn','',...
            'WindowButtonDownFcn','',...
            'WindowButtonUpFcn','');
end

% Update structure
guidata(hObject,handles);


function edit_start_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start as text
%        str2double(get(hObject,'String')) returns contents of edit_start as a double

try
    % Get the start and stop time in duration values
    split = strsplit(handles.edit_start.String, ':');
    temp_start = duration(str2double(split));
    
    stop = handles.temp.stop.XData(1);
    
    if ~isnan(temp_start) && ~sum(contains(split,[",","."]))
        % Check if it is bigger than the stop time
        if temp_start > stop
            % Adjust the edit box
            set(handles.edit_start,'string',char(stop))
            
            % Adjust the start time
            start = stop;
            
        else
            start = temp_start;
        end
        
        % Adjust the start line
        set(handles.temp.start,'XData',[start start],...
            'YData',ylim(handles.axes_ECG),...
            'visible','on');
        
        % Adjust the patch
        Ylim = ylim(handles.axes_ECG);
        
        handles.temp.patch.XData([1 4]) = seconds([start start]);
        set(handles.temp.patch,'YData',[Ylim(1) Ylim(1) Ylim(2) Ylim(2)],...
            'visible','on')
        
        % Adjust the start value and edit box
        handles.data.start_analysis = start;
        set(handles.edit_start,'string',char(start))
        
        % Adjust the duration value and edit box
        handles.data.duration_analysis = stop-start;
        set(handles.edit_duration,'string',char(stop-start))
    else
        % Set start to previous value
        set(handles.edit_start,'String',char(handles.data.start_analysis));
        
        % Throw an error dialog
        errordlg('Start time should be in the "hh:mm:ss" format.','Format error')
    end
catch
    % Set start to previous value
    set(handles.edit_start,'String',char(handles.data.start_analysis));
    
    % Throw an error dialog
    errordlg('Start time should be in the "hh:mm:ss" format.','Format error')
end

% Update structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_start_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_duration_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double

try
    % Get the start and duration in duration values
    start = handles.temp.start.XData(1);
    
    split = strsplit(handles.edit_duration.String, ':');
    temp_duration = duration(str2double(split));
    
    if ~isnan(temp_duration) && ~sum(contains(split,[",","."]))
        % Compute the new stop
        temp_stop = start + temp_duration;
        
        % Check if it is beyond the length of the signal
        if temp_stop > handles.data.duration_recording
            % Define stop
            stop = handles.data.duration_recording;
            
            % Adjust the duration value and edit box
            handles.data.duration_analysis = stop-start;
            set(handles.edit_duration,'string',char(stop - start))
            
        else
            stop = temp_stop;
            
            % Adjust the duration value and edit box
            handles.data.duration_analysis = temp_duration;
            set(handles.edit_duration,'string',char(temp_duration))
        end
        
        % Adjust the patch
        handles.temp.patch.XData([2 3]) = seconds([stop stop]);
        
        % Adjust the stop line
        handles.temp.stop.XData = [stop stop];
        
    else
        % Set duration to previous value
        set(handles.edit_duration,'String',char(handles.data.duration_analysis));
        
        % Throw an error dialog
        errordlg('Duration should be in the "hh:mm:ss" format.','Format error')
    end
catch
    % Set duration to previous value
    set(handles.edit_duration,'String',char(handles.data.duration_analysis));
    
    % Throw an error dialog
    errordlg('Duration should be in the "hh:mm:ss" format.','Format error')
end

% Update structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_detect_peaks.
function pushbutton_detect_peaks_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_detect_peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')

% Reset the manual corrections
set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)

% Execute the radiobutton callback
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Define sampling frequency
fs = handles.data.fs;

% Define the entire time
temp_time = seconds((1:size(handles.data.signal.filtered,1))/fs);
temp_time.Format = 'hh:mm:ss';

% Define start and stop
temp_start = handles.data.start_analysis;
temp_stop = temp_start + handles.data.duration_analysis;

% Define the to be analysed signal
start = find(temp_time >= temp_start,1,'first');
stop = find(temp_time < temp_stop,1,'last');

signal = handles.data.signal.filtered(start:stop,:);

% Define the time selection
time = seconds((start:stop)/fs);

% Detect the R-peaks
[R_peak, ~, parameters, check] = peak_detection(handles.data.parameters.Rpeak,signal,fs,1);

if check
    % Adjust the R-peaks
    handles.data.R = [];
    for ii = 1:handles.data.channels
        handles.data.R{1,ii} = time(R_peak{ii});
    end
    
    % Adjust the parameters
    handles.data.parameters.Rpeak = parameters;
    
    % Adjust the plots
    handles = update_plots(handles);
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
    
    % Enable add, adjust and delete
    set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','on')
    
    % Enable the add context menu
    set(handles.context_add,'enable','on')
    
    % Enable the menu's
    set(handles.menu_export_results,'enable','on')
    set(handles.menu_save_results,'enable','on')
end

% Change the pointer
set(handles.figure_main,'Pointer','arrow')

% Update structure
guidata(hObject,handles);


% --- Executes on button press in pushbutton_load_peaks.
function pushbutton_load_peaks_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_load_peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Reset the manual corrections
set(findall(handles.uipanel_R_peak_correction,'Style','radiobutton'),'value',0)

% Execute the radiobutton callback
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Select where you want to load from
answer = questdlg('Where do you want to load the peaks from?', ...
    'Load Peaks', ...
    'From File','From Workspace','From File');

% Set check
check = 1;
try
    switch answer
        case 'From File'
            % Select the file
            [FileName,PathName] = uigetfile({'*.mat'},'Select the R-peaks');
            
            % If the window was not closed, check the file format
            if ~isnumeric(FileName) && ~isnumeric(PathName)
                
                % Get the complete file
                file = strcat(PathName,FileName);
                
                % Import the file
                data = (importdata(file));
                
            end
        case 'From Workspace'
            % Select the variables from the workspace
            vars = evalin('base','who');
            
            % Make a selection if multiple variables are present
            if length(vars) > 1
                [idx,~] = listdlg('PromptString','Select the R-peak locations:',...
                    'ListString',vars,...
                    'SelectionMode','single');
            else
                idx = 1;
            end
            
            % Load the variable
            data = evalin('base',vars{idx});
            
        otherwise
            check = 0;
    end
catch
    check = 0;
end

if check == 1
    try
        if isa(data,'struct') % If var is a structure, loop through the structure
            % Select the correct file
            [R_peak_temp,~] = loop_through_structure(data);
        else % If var is a variable, select the variable
            R_peak_temp = data;
        end
        
        if isnumeric(R_peak_temp)
            
            % Set the signal in the correct form (collumn vectors)
            [R_temp, K_temp] = size(R_peak_temp);
            if K_temp < R_temp
                R_peak_temp = R_peak_temp';
                R = K_temp;
            else
                R = R_temp;
            end
            
            % Ask in which units the R-peaks are
            [idx,check] = listdlg('PromptString','Select the R-peak units:',...
                'ListString',{'Samples','ms','s'},...
                'SelectionMode','single');
            
            if check
                if R > 1
                    % Pre-allocate
                    R_peak = cell(1,R);
                    
                    % Loop over the leads
                    for ii = 1:R
                        if idx == 1
                            R_peak{ii} = R_peak_temp(ii,:);
                        elseif idx == 2
                            R_peak{ii} = round((R_peak_temp(ii,:)/1000)*handles.data.fs);
                        elseif idx == 3
                            R_peak{ii} = round(R_peak_temp(ii,:)*handles.data.fs);
                        end
                    end
                else
                    if handles.data.channels == 1
                        R_peak{1} = R_peak_temp;
                    else
                        % Pre-allocate
                        R_peak = cell(1,handles.data.channels);
                        
                        % Loop over the leads
                        for ii = 1:handles.data.channels
                            R_peak{ii} = R_peak_temp;
                        end
                    end
                end
            else
                return
            end
        elseif iscell(R_peak_temp)
            % Store the cell
            R_peak = R_peak_temp;
            
            % Get the size of the cell
            %             R = length(R_peak);
        elseif isduration(R_peak_temp)
            % Store as cell
            R_peak = {R_peak_temp};
        end
        
        % Check if the R-peak locations are in ms, s or samples
        if ~isduration(R_peak{1})
            
            % Define sampling frequency
            fs = handles.data.fs;
            
            % Define the entire time
            temp_time = seconds((1:size(handles.data.signal.filtered,1))/fs);
            temp_time.Format = 'hh:mm:ss';
            
            % Define start and stop
            temp_start = handles.data.start_analysis;
            temp_stop = temp_start + handles.data.duration_analysis;
            
            % Define the to be analysed signal
            start = find(temp_time >= temp_start,1,'first');
            stop = find(temp_time < temp_stop,1,'last');
            
            % Define the time selection
            time = seconds((start:stop)/fs);
            
            % Adjust the R-peaks
            handles.data.R = [];
            for ii = 1:handles.data.channels
                handles.data.R{1,ii} = time(R_peak{ii});
            end
        else % Isduration
            % Adjust the R-peaks
            handles.data.R = [];
            for ii = 1:handles.data.channels
                handles.data.R{1,ii} = R_peak{ii};
            end
        end
        
        % Adjust the plots
        handles = update_plots(handles);
        
        % Adjust the y-limits
        handles = adjust_y_limits(handles);
        
        % Enable add, adjust and delete
        set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','on')
        
        % Enable the add context menu
        set(handles.context_add,'enable','on')
        
        % Enable the menu's
        set(handles.menu_export_results,'enable','on')
        set(handles.menu_save_results,'enable','on')
    catch
    end
end

% Change the pointer
set(handles.figure_main,'Pointer','arrow')
drawnow;

% Update structure
guidata(hObject,handles);


% --- Executes on button press in pushbutton_cross_correlation.
function pushbutton_cross_correlation_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_cross_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Define the looping variable
if listbox_value == 0
    loop_var = 1:handles.data.channels;
else % Only the selection
    loop_var = listbox_value;
end

% Define sampling frequency
fs = handles.data.fs;

% Define the signal
signal = handles.data.signal.filtered;

% Define the time
time = seconds((1:size(signal,1))/fs);
time.Format = 'hh:mm:ss';

% Pre-allocate
window = 0.15;

% Loop over all the R-peaks
for ii = loop_var
    % Define a temporary R-peak variable
    R_peak_temp = round(fs*seconds(handles.graph.R.(handles.labels{ii}).XData));
    
    % Avoid trouble
    while R_peak_temp(1)-round(window*fs) <= 0
        R_peak_temp(1) = [];
    end
    while R_peak_temp(end)+round(window*fs) > length(signal)
        R_peak_temp(end) = [];
    end
    
    % Define the size of the heartbeat variable
    R = round(window*fs)*2+1;
    K = length(R_peak_temp);
    
    % Select all heartbeats
    QRS = zeros(R,K);
    for iii = 1:K
        % Slide 200ms forward and backward
        QRS(:,iii) = signal(R_peak_temp(iii)-round(window*fs):R_peak_temp(iii)+round(window*fs),ii);
    end
    
    % Normalize the beats
    QRS = zscore(QRS);
    
    % Get the ammount of positive R-peaks
    nr_pos = sum(sign(signal(R_peak_temp,ii)));
    
    % Predefine the shift
    R_shift = 0;
    
    try
        % Compute the average positive R-beat
        QRS_avg_pos = trimmean(QRS(:,sign(signal(R_peak_temp,ii))==1)',30)';
        
        % Compute the average "negative" R-beat
        QRS_avg_neg = trimmean(QRS(:,sign(signal(R_peak_temp,ii))==-1)',30)';
        
        % Select the positive or negative R-peak
        [selection, R_shift] = QRS_selection(QRS_avg_pos,QRS_avg_neg,window,nr_pos,fs);
        
    catch
        if nr_pos > 0
            selection = 1;
        else
            selection = -1;
        end
    end
    
    % Compute the average heartbeat
    QRS_avg = trimmean(QRS(:,sign(signal(R_peak_temp,ii)) == selection)',30)';
    
    % Pre-allocate
    for iii = 1:K
        % Compute cross1-correlation
        [QRS_corr,lag] = xcorr(QRS(:,iii), QRS_avg);
        
        % Adjust the R-peak position
        [~,I] = max(QRS_corr);
        R_peak_temp(iii) = R_peak_temp(iii)+lag(I)+R_shift;
    end
    
    % Store the R-peaks and compute RR-intervals
    set(handles.graph.R.(handles.labels{ii}),'XData',time(round(unique(R_peak_temp))),...
        'YData',signal(round(unique(R_peak_temp)),ii))
end

% Check if the manual corrections are selected and update the R-peaks model
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Adjust the plots
handles = update_plots(handles);

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Enable the save menu's
set(handles.menu_export_results,'enable','on')
set(handles.menu_save_results,'enable','on')

% Change the pointer
set(handles.figure_main,'Pointer','arrow')
drawnow;

% Update structure
guidata(hObject,handles)


% --- Executes on button press in pushbutton_peak_conversion.
function pushbutton_peak_conversion_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_peak_conversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Define the looping variable
if listbox_value == 0
    loop_var = 1:handles.data.channels;
else % Only the selection
    loop_var = listbox_value;
end

% Define sampling frequency
fs = handles.data.fs;

% Define the signal
signal = handles.data.signal.filtered;

% Define the time
time = seconds((1:size(signal,1))/fs);
time.Format = 'hh:mm:ss';

% Loop over all the R-peaks
for ii = loop_var
    % Define a temporary R-peak variable
    R_peak_temp = round(fs*seconds(handles.graph.R.(handles.labels{ii}).XData));
    
    % Converge...
    for iii = 1:length(R_peak_temp)
        % ... forward
        while abs(signal(R_peak_temp(iii),ii)) > abs(signal(R_peak_temp(iii)-1,ii)) && ...
                abs(signal(R_peak_temp(iii),ii)) < abs(signal(R_peak_temp(iii)+1,ii))
            R_peak_temp(iii) = R_peak_temp(iii)+1;
        end
        
        % ... backward
        while abs(signal(R_peak_temp(iii),ii)) < abs(signal(R_peak_temp(iii)-1,ii)) && ...
                abs(signal(R_peak_temp(iii),ii)) > abs(signal(R_peak_temp(iii)+1,ii))
            R_peak_temp(iii) = R_peak_temp(iii)-1;
        end
    end
    
    % Store the R-peaks and compute RR-intervals
    set(handles.graph.R.(handles.labels{ii}),'XData',time(round(unique(R_peak_temp))),...
        'YData',signal(round(unique(R_peak_temp)),ii))
end

% Check if the manual corrections are selected and update the R-peaks model
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Adjust the plots
handles = update_plots(handles);

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Change the pointer
set(handles.figure_main,'Pointer','arrow')
drawnow;

% Update structure
guidata(hObject,handles)


% --- Executes on button press in pushbutton_maximum_search.
function pushbutton_maximum_search_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_maximum_search (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the pointer
set(handles.figure_main,'Pointer','watch')
drawnow;

% Ask if they want a max, min or abs max search
choice = questdlg('Select the type of search.','Maximum search','Max','Min','Abs Max','Max');

if ~strcmp(choice,'')
    % Get the listbox value
    listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...
    
    % Define the looping variable
    if listbox_value == 0
        loop_var = 1:handles.data.channels;
    else % Only the selection
        loop_var = listbox_value;
    end
    
    % Define sampling frequency
    fs = handles.data.fs;
    
    % Define the signal
    signal = handles.data.signal.filtered;
    
    % Define the time
    time = seconds((1:size(signal,1))/fs);
    time.Format = 'hh:mm:ss';
    
    % Define the window
    window = round(0.1*fs);
    
    % Loop over all the R-peaks
    for ii = loop_var
        try
            % Define a temporary R-peak variable
            R_peak_temp = round(fs*seconds(handles.graph.R.(handles.labels{ii}).XData));
            
            % Find maximum
            for iii = 1:length(R_peak_temp)
                % Get the samples from the R-peak minus and plus the window
                p = signal(max(1,R_peak_temp(iii)-window):min(R_peak_temp(iii)+window,length(signal)),ii);
                
                % Do the same but for the time locations
                t = max(1,R_peak_temp(iii)-window):min(R_peak_temp(iii)+window,length(signal));
                
                % Do a maximum search based on the selection
                switch choice
                    case 'Max'
                        [~,ip] = max(p);
                    case 'Min'
                        [~,ip] = min(p);
                    case 'Abs Max'
                        [~,ip] = max(abs(p));
                end
                
                % Adjust the R-peak
                R_peak_temp(iii) = t(ip(end));
            end
            
            % Store the R-peaks and compute RR-intervals
            set(handles.graph.R.(handles.labels{ii}),'XData',time(round(unique(R_peak_temp))),...
                'YData',signal(round(unique(R_peak_temp)),ii))
        catch
        end
    end
    
    % Check if the manual corrections are selected and update the R-peaks model
    handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
    
    % Adjust the plots
    handles = update_plots(handles);
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end

% Change the pointer
set(handles.figure_main,'Pointer','arrow')
drawnow;

% Update structure
guidata(hObject,handles)


% --- Executes on button press in radiobutton_add, radiobutton_adjust, radiobutton_delete.
function handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_add, radiobutton_adjust and radiobutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the toggle tools
% (Keeping this on will give problems when changing the WindowButtonDownFcn)
% set(findall(handles.figure_main, 'Type','uitoggleTool'), 'enable', 'off')
zoom off

% Check if one of the radiobuttons is one
if get(handles.radiobutton_add,'Value') || get(handles.radiobutton_adjust,'Value') || get(handles.radiobutton_delete,'Value')
    % Remove all window functions
    set(handles.figure_main,'WindowButtonMotionFcn','',...
        'WindowButtonDownFcn','',...
        'WindowButtonUpFcn','');
    
    % Reset the other radiobuttons
    % This will only work if one of the radiobuttons is selected
    try
        switch hObject.String
            case 'Add'
                set([handles.radiobutton_adjust handles.radiobutton_delete],'value',0)
            case 'Adjust'
                set([handles.radiobutton_add handles.radiobutton_delete],'value',0)
            case 'Delete'
                set([handles.radiobutton_add handles.radiobutton_adjust],'value',0)
        end
    catch
    end
    
    if get(handles.radiobutton_add,'Value')
        
        % Create a cursor line
        handles.temp.cursor_line = line(handles.axes_ECG,seconds([NaN NaN]), ylim(handles.axes_ECG), ...
            'Color', 'black',...
            'Parent', handles.axes_ECG);
        
        % Create a temporary R-peak for every lead
        for ii = 1:handles.data.channels
            handles.temp.R(ii) = plot(handles.axes_ECG,seconds(NaN),handles.data.signal.filtered(1),...
                'og',...
                'markerfacecolor','g',...
                'linewidth',3);
        end
        
        % Set the window buttonmotion and down function
        set(handles.figure_main,'WindowButtonMotionFcn',{@drag_new_R_peak,handles},...
            'WindowButtonDownFcn',{@define_new_R_peak,handles});
        
    elseif get(handles.radiobutton_adjust,'Value')
        
        % Create a temporary R-peak
        for ii = 1:handles.data.channels
            handles.temp.R(ii) = plot(handles.axes_ECG,seconds(NaN),handles.data.signal.filtered(1),...
                'oc',...
                'markerfacecolor','c',...
                'linewidth',3);
        end
        
        % Set up cursor line
        handles.temp.cursor_line = line(handles.axes_ECG,seconds([NaN NaN]), ylim(handles.axes_ECG), ...
            'Color', 'black',...
            'Parent', handles.axes_ECG);
        
        % Set the window buttondown and up function
        set(handles.figure_main,'WindowButtonMotionFcn',{@select_R_peak,handles},...
            'WindowButtonDownFcn',{@define_R_peak,handles},...
            'WindowButtonUpFcn',{@let_R_peak_go,handles});
        
    elseif get(handles.radiobutton_delete,'Value')
        
        % Create a temporary R-peak for every lead
        for ii = 1:handles.data.channels
            handles.temp.R(ii) = plot(handles.axes_ECG,seconds(NaN),handles.data.signal.filtered(1),...
                'or',...
                'markerfacecolor','r',...
                'linewidth',3);
        end
        
        % Set the window buttonmotion and down function
        set(handles.figure_main,'WindowButtonMotionFcn',{@select_R_peak,handles},...
            'WindowButtonDownFcn',{@change_buttonmotionfcn,handles},...
            'WindowButtonUpFcn',{@delete_R_peak,handles});
    end
    
else
    
    % Set the window buttonmotion and down function
    set(handles.figure_main,'WindowButtonMotionFcn','',...
        'WindowButtonDownFcn','',...
        'WindowButtonUpFcn','');
    
    % Adjust the model
    if ~isempty(handles.data.R)
        for ii = 1:handles.data.channels
            % Find the R-peaks in the model in the search interval
            temp = find(handles.data.R{ii} >= handles.data.start_analysis);
            idx = temp(handles.data.R{ii}(temp) <= handles.data.start_analysis + handles.data.duration_analysis);
            
            if length(idx) == length(handles.graph.R.(handles.labels{ii}).XData)
                handles.data.R{ii}(idx) = handles.graph.R.(handles.labels{ii}).XData;
            elseif handles.graph.R.(handles.labels{ii}).XData(1) < handles.data.R{ii}(1)
                handles.data.R{ii} = [handles.graph.R.(handles.labels{ii}).XData handles.data.R{ii}];
            elseif handles.graph.R.(handles.labels{ii}).XData(1) > handles.data.R{ii}(end)
                handles.data.R{ii} = [handles.data.R{ii} handles.graph.R.(handles.labels{ii}).XData];
            else
                handles.data.R{ii} = [handles.data.R{ii}(1:idx(1)-1) handles.graph.R.(handles.labels{ii}).XData handles.data.R{ii}(idx(end)+1:end)];
            end
        end
    end
    
    % Delete temporary graphics
    try
        delete(handles.temp.R)
        delete(handles.temp.cursor_line)
    catch
    end
end

% Update structure
guidata(hObject,handles)


function handles = edit_range_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_range as text
%        str2double(get(hObject,'String')) returns contents of edit_range as a double

try
    % Get the range in duration values
    split = strsplit(handles.edit_range.String, ':');
    range = duration(str2double(split));
    
    if ~isnan(range) && ~sum(contains(split,[",","."]))
        
        % Get the start of the analysis period in duration values
        an_start = handles.data.start_analysis;
        
        % Get the width of the analysis window in duration values
        an_width = handles.data.duration_analysis;
        
        % Get the stop of the analysis period in duration values
        an_stop = an_start + an_width;
        
        if range >= an_width
            % Change the edit box
            set(handles.edit_range,'string',char(an_width))
            
            % Adjust the slider
            set(handles.slider_range,'enable','off',...
                'max',round(handles.data.fs*seconds(an_width)),...
                'value',round(handles.data.fs*seconds(an_width)))
            
            % Adjust range
            range = an_width;
            
        elseif range < seconds(1)
            % Set the string
            set(handles.edit_range,'string','00:00:01')
            
            % Adjust range
            range = seconds(1);
            
        else
            % Set the string
            set(handles.edit_range,'string',char(range));
            
            % Get a new maximum for the slider
            new_max = round(handles.data.fs*seconds(an_width-range));
            
            % Adjust the value if it is outside the new range
            if handles.slider_range.Value > new_max
                set(handles.slider_range,'Value',new_max)
            end
            
            % Get the new step size
            n = round(handles.data.fs*seconds(range))/new_max;
            
            % Adjust the slider
            set(handles.slider_range,'max',new_max,...
                'sliderstep',[min([1 n]) min([1 2*n])],...
                'enable','on')
        end
        
        % Get the slider value
        slider_value = get(handles.slider_range,'Value');
        
        % Adjust the x-limits
        stop = min([an_start+seconds(slider_value/handles.data.fs)+range an_stop]);
        
        if stop == an_stop
            if strcmp(handles.slider_range.Enable,'on')
                start = stop-range;
            else
                start = an_start;
            end
        else
            start = an_start+seconds(slider_value/handles.data.fs);
        end
        
        % Adjust the x-limits
        xlim(handles.axes_ECG,[start stop])
        
        % Adjust the y-limits
        handles = adjust_y_limits(handles);
        
        % Set the model range value
        handles.data.range = duration(str2double(strsplit(handles.edit_range.String, ':')));
    else
        % Reset the string
        set(handles.edit_range,'string',char(handles.data.range));
        
        % Throw an error dialog
        errordlg('Range should be in the "hh:mm:ss" format.','Format error')
    end
catch
    % Reset the string
    set(handles.edit_range,'string',char(handles.data.range));
    
    % Throw an error dialog
    errordlg('Range should be in the "hh:mm:ss" format.','Format error')
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_range_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plus.
function pushbutton_plus_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the interval width in duration values
temp = duration(str2double(strsplit(handles.edit_range.String, ':')));

% Add one second
temp = temp+seconds(1);

% Update range string
set(handles.edit_range,'String',char(temp))

% Update range
handles = edit_range_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_plus.
function pushbutton_plus_ButtonDownFcn(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get clicktype
clickType = get(handles.figure_main, 'SelectionType');

if strcmp(handles.pushbutton_plus.Enable,'on') && strcmp(clickType, 'alt')
    % Get the start of the analysis period in duration values
    start = handles.data.start_analysis;
    
    % Get the width of the analysis window in duration values
    an_width = handles.data.duration_analysis;
    
    stop = start+an_width;
    
    % Set the range edit box
    set(handles.edit_range,'string',char(an_width))
    
    % Disable the slider, (re)set the slider maximum and change
    % the value
    set(handles.slider_range,'enable','off',...
        'max',round(handles.data.fs*seconds(an_width)),...
        'value',0)
    
    % Adjust the x-limits
    xlim(handles.axes_ECG,[start stop])
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
    
    % Set the model
    handles.data.range = an_width;
    
    try
        if strcmp(handles.temp.patch.Visible,'on')
            Ylim = handles.axes_ECG.YLim;
            set(handles.temp.patch,'YData',[Ylim(1) Ylim(1) Ylim(2) Ylim(2)])
            set(handles.temp.start,'YData',Ylim)
            set(handles.temp.stop,'YData',Ylim)
        end
    catch
    end
    
    % Reset the zoom
    zoom 'reset'
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_minus.
function pushbutton_minus_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the interval width in duration values
temp = duration(str2double(strsplit(handles.edit_range.String, ':')));

% Add one second
temp = temp-seconds(1);

% Update range string
set(handles.edit_range,'String',char(temp))

% Update range
handles = edit_range_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_range_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to slider_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_range_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to slider_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listbox_channels.
function listbox_channels_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to listbox_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_channels

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Get the selection type
selectiontype = get(handles.figure_main,'SelectionType');

% If double click
if strcmp(selectiontype,'open') && listbox_value > 0
    % Get all names
    all_names = handles.listbox_channels.String;
    
    % Get the selection in the listbox
    idx = listbox_value + 1;
    
    % Get the old name
    L = length(all_names{idx});
    string = all_names{idx};
    old_name = string(find(string=='"',1,'last')+2:L-14);
    
    % Get the new name
    new_name = rename_channel('channel_name',old_name);
    
    if ~isempty(new_name)
        % Adjust the listbox
        color = handles.color(listbox_value,:);
        colorStr = sprintf('%d,',int16(255*color));
        all_names{idx} = ['<HTML><FONT color="rgb(' colorStr ')">' new_name '</Font></html>'];
        set(handles.listbox_channels,'string',all_names)
    end
else
    % Define the looping variable
    if listbox_value == 0
        loop_var = 1:handles.data.channels;
    else % Only the selection
        loop_var = listbox_value;
    end
    
    % Get the original signal checkbox value
    original = get(handles.checkbox_show_original,'value');
    
    % Make all signals invisible
    set(handles.axes_ECG.Children,'Visible','off')
    set(handles.axes_tachogram.Children,'Visible','off')
    
    % Make the selected signals visible
    for ii = loop_var
        % Filtered ECG signals
        set(handles.graph.signal.filtered.(handles.labels{ii}),'visible','on')
        
        % Original ECG signals
        if original
            set(handles.graph.signal.original.(handles.labels{ii}),'visible','on')
        end
        
        % R-peaks + tachogram
        if ~isempty(handles.data.R)
            set(handles.graph.R.(handles.labels{ii}),'visible','on')
            set(handles.graph.RR.(handles.labels{ii}),'visible','on')
        end
    end
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end

% Update structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox_channels_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to listbox_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_units_tachogram.
function popupmenu_units_tachogram_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to popupmenu_units_tachogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_units_tachogram contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_units_tachogram

% Get the contents
contents = cellstr(get(handles.popupmenu_units_tachogram,'String'));
value = get(handles.popupmenu_units_tachogram,'Value');
units = contents{value};

% Adjust the model
handles.data.units = value;

% Change the label
ylabel(handles.axes_tachogram, units)

% Define the RR-intervals
if ~isempty(handles.data.R)
    
    % Loop over the channels
    for ii = 1:handles.data.channels
        % Compute the RR-interval
        RR_int = seconds(diff(handles.graph.R.(handles.labels{ii}).XData))*1000;
        
        switch units
            case 'HR (bpm)'
                % Change the data
                RR_int = 60000./RR_int;
        end
        
        % Adjust the tachogram
        handles.graph.RR.(handles.labels{ii}).YData = RR_int;
    end
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end

% Update structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_units_tachogram_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to popupmenu_units_tachogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_fix_ylimits.
function checkbox_fix_ylimits_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to checkbox_fix_ylimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fix_ylimits

% Get the checkbox value
checkbox_value = get(handles.checkbox_fix_ylimits,'Value');

if ~checkbox_value
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end

% Update structure
guidata(hObject,handles);


%% Extra functions %%
function Slide(hObject,eventdata) %#ok

% Get the handles
handles = guidata(eventdata.AffectedObject);

if strcmp(handles.slider_range.String,'normal')
    % Get the range in duration values
    range = duration(str2double(strsplit(handles.edit_range.String, ':')));
    
    % Get the start of the analysis period in duration values
    an_start = handles.data.start_analysis;
    
    % Get the width of the analysis window in duration values
    an_width = handles.data.duration_analysis;
    
    % Get the stop of the analysis period in duration values
    an_stop = an_start + an_width;
    
    % Get the slider value
    slider_value = get(eventdata.AffectedObject,'Value');
    
    % Adjust the x-limits
    stop = min([an_start+seconds(slider_value/handles.data.fs)+range an_stop]);
    
    if stop == an_stop
        if strcmp(handles.slider_range.Enable,'on')
            start = stop-range;
        else
            start = an_start;
        end
    else
        start = an_start+seconds(slider_value/handles.data.fs);
    end
    
    % Adjust the x-limits
    xlim(handles.axes_ECG,[start stop])
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
    
    % Draw everything now. This makes sure that the function will execute, also
    % when the zoom is enabled
    if strcmp(handles.zoom.Enable,'on')
        drawnow
    end
end


%% Helper functions analysis period

function drag_start_line(hObject,eventdata,handles) %#ok
% windowbuttonmotionfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the candidate position
candidate_pos = mousepos(1);

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) && ...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) && ...
        mousepos(1,2) <= Ylim(2)
    
    % Adjust the mouseposition
    if candidate_pos <= 0 % Past the left limit
        new_pos = 0;
    elseif candidate_pos >= seconds(handles.temp.stop.XData(1))-1 % Past the right limit
        new_pos = seconds(handles.temp.stop.XData(1))-1;
    else % Perfect
        new_pos = candidate_pos;
    end
    
    % Set the new start time
    start = seconds(new_pos);
    start.Format = 'hh:mm:ss';
    set(handles.edit_start,'string',char(start));
    
    % Make the starting line visible
    set(handles.temp.start,'XData',[start start],...
        'YData',Ylim,...
        'visible','on')
    
    % Adjust the patch
    handles.temp.patch.XData([1 4]) = [new_pos new_pos];
    handles.temp.patch.YData = [Ylim(1) Ylim(1) Ylim(2) Ylim(2)];
end

% Adjust the duration edit box
start = duration(str2double(strsplit(handles.edit_start.String, ':')));
stop = handles.temp.stop.XData(1);
set(handles.edit_duration,'string',char(stop-start))

% Update structure
guidata(hObject,handles);


function define_start_line(hObject,eventdata,handles) %#ok
% windowbuttondownfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the candidate position
candidate_pos = mousepos(1);

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    % Get the x-location
    new_pos = candidate_pos;
    
    % Make the patch visible
    set(handles.temp.patch,'XData',[new_pos new_pos new_pos new_pos],...
        'YData',[Ylim(1) Ylim(1) Ylim(2) Ylim(2)],...
        'visible','on')
    
    % Make the stopping line visible
    set(handles.temp.stop,'XData',[new_pos new_pos],...
        'YData',Ylim,...
        'visible','on')
    
    % Set the window button motion function
    set(handles.figure_main,'WindowButtonMotionFcn',{@drag_stop_line,handles});
end

% Update structure
guidata(hObject,handles);


function drag_stop_line(hObject,eventdata,handles) %#ok
% windowbuttonmotionfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the candidate position
candidate_pos = mousepos(1);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check the constraints
if candidate_pos <= seconds(handles.temp.start.XData(1)) % Past the start line
    new_pos = seconds(handles.temp.start.XData(1))+1;
elseif candidate_pos >= seconds(handles.data.duration_recording) % Past the end of the recording
    new_pos = seconds(handles.data.duration_recording);
else % Perfect
    new_pos = mousepos(1);
end

% Get the location of the start line
temp_start = duration(str2double(strsplit(handles.edit_start.String, ':')));
temp_start.Format = 'hh:mm:ss';

% Get the location of the stop line
temp_stop = seconds(new_pos);
temp_stop.Format = 'hh:mm:ss';

% Adjust the start edit box
set(handles.edit_start,'string',char(temp_start));

% Adjust the duration edit box
set(handles.edit_duration,'string',char(temp_stop-temp_start))

% Adjust the stop line
set(handles.temp.stop,'XData',[temp_stop temp_stop],...
    'YData',Ylim);

% Adjust the patch
handles.temp.patch.XData([2 3]) = [new_pos new_pos];
handles.temp.patch.YData = [Ylim(1) Ylim(1) Ylim(2) Ylim(2)];

% Update structure
guidata(hObject,handles);


function define_stop_line(hObject,eventdata,handles) %#ok
% windowbuttonupfcn

switch char(handles.figure_main.WindowButtonDownFcn{1,1})
    case 'define_start_line'
        
        if strcmp(char(handles.figure_main.WindowButtonMotionFcn{1,1}),'drag_stop_line')
            
            % Remove the window button motion function
            set(handles.figure_main,'WindowButtonMotionFcn','');
            
            % Enable the analysis period edit boxes
            set(findall(handles.uipanel_analysis_period,'style','edit'),'enable','on')
            
            % Enable the pointer manager
            iptPointerManager(handles.figure_main,'enable');
            
            % Change the pointer behavior for the lines
            pointerBehavior.enterFcn    = [];
            pointerBehavior.exitFcn     = [];
            pointerBehavior.traverseFcn = @overline;
            
            iptSetPointerBehavior(handles.temp.start, pointerBehavior);
            iptSetPointerBehavior(handles.temp.stop, pointerBehavior);
            
            % Change the pointer behavior for the patch
            pointerBehavior.traverseFcn = @cross1;
            
            iptSetPointerBehavior(handles.temp.patch, pointerBehavior);
            
            % Define new windowbuttondown function
            set(handles.figure_main,'WindowButtonDownFcn',{@select_line,handles});
        end
    case 'select_line'
        % Remove the window button motion function
        set(handles.figure_main,'WindowButtonMotionFcn','');
        
        % Enable the pointer manager
        iptPointerManager(handles.figure_main,'enable');
        
        % Change the pointer behavior for the lines
        pointerBehavior.enterFcn    = [];
        pointerBehavior.exitFcn     = [];
        pointerBehavior.traverseFcn = @overline;
        
        iptSetPointerBehavior(handles.temp.start, pointerBehavior);
        iptSetPointerBehavior(handles.temp.stop, pointerBehavior);
        
        % Change the pointer behavior for the patch
        pointerBehavior.traverseFcn = @cross1;
        
        iptSetPointerBehavior(handles.temp.patch, pointerBehavior);
end

% Update structure
guidata(hObject,handles);


function drag_patch(hObject,eventdata,handles) %#ok
% windowbuttonmotionfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-position
x_pos = mousepos(1);

% Get the current position
current_pos = [ seconds(handles.temp.start.XData(1)) seconds(handles.temp.stop.XData(1)) ];

% Get the candidate position
candidate_pos = [x_pos-handles.temp.start_to_mouse x_pos+handles.temp.stop_to_mouse];

% Check the constraints
if candidate_pos(1) < 0
    new_pos(1) = 0;
    new_pos(2) = diff(current_pos);
    
elseif candidate_pos(2) > seconds(handles.data.duration_recording)
    new_pos(2) = seconds(handles.data.duration_recording);
    new_pos(1) = new_pos(2) - diff(current_pos);
else
    new_pos = candidate_pos;
end

% Set the new start time
start = seconds(new_pos(1));
start.Format = 'hh:mm:ss';
set(handles.edit_start,'string',char(start));

% Adjust the start line
set(handles.temp.start,'XData',[start start])

% Adjust the patch
handles.temp.patch.XData([1 4]) = [new_pos(1) new_pos(1)];

% Get the location of the stop line
stop = seconds(new_pos(2));
stop.Format = 'hh:mm:ss';

% Adjust the stop line
set(handles.temp.stop,'XData',[stop stop]);

% Adjust the patch
handles.temp.patch.XData([2 3]) = [new_pos(2) new_pos(2)];

% Update structure
guidata(hObject,handles);


function select_line(hObject,eventdata,handles) %#ok
% windowbuttondownfcn

% Check if the pointer is a double arrow
if strcmp(handles.figure_main.Pointer,'left')
    % Get the clicked object
    co = gco;
    
    % Get the current objects x-location
    idx = round(seconds(co.XData(1)),2);
    
    % Disable the pointer manager
    iptPointerManager(handles.figure_main,'disable');
    
    % Fix it to left
    set(handles.figure_main,'Pointer','left')
    
    % Check which line is selected
    if idx == round(handles.temp.patch.XData(1),2)
        % Set the window button motion function
        set(handles.figure_main,'WindowButtonMotionFcn',{@drag_start_line,handles});
    elseif idx == round(handles.temp.patch.XData(2),2)
        % Set the window button motion function
        set(handles.figure_main,'WindowButtonMotionFcn',{@drag_stop_line,handles});
    end
elseif strcmp(handles.figure_main.Pointer,'fleur')
    % Disable the pointer manager
    iptPointerManager(handles.figure_main,'disable');
    
    % Fix it to fleur
    set(handles.figure_main,'Pointer','fleur')
    
    % Get the mouse position
    temp = get(handles.axes_ECG,'CurrentPoint');
    
    % Distance of the start line to the mouse position
    handles.temp.start_to_mouse = temp(1) - seconds(handles.temp.start.XData(1));
    
    % Distance of the stop line to the mouse position
    handles.temp.stop_to_mouse = seconds(handles.temp.stop.XData(1)) - temp(1);
    
    % Set the window button motion function
    set(handles.figure_main,'WindowButtonMotionFcn',{@drag_patch,handles});
end

% Update structure
guidata(hObject,handles);


%% Helper functions add, delete and adjust radiobuttons

function select_R_peak(hObject,eventdata,handles) %#ok
% windowbuttonmotionfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));
Xwidth = Xlim(2)-Xlim(1);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2)-Ylim(1);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    if ~listbox_value % If all leads are displayed
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    % If the delete radiobutton is selected or the adjust button and one
    % lead
    if handles.radiobutton_delete.Value
        for ii = loop_var
            % Get the closest R-peak
            [~,idx] = min(((seconds(handles.graph.R.(handles.labels{ii}).XData)-mousepos(1))./Xwidth).^2+...
                ((handles.graph.R.(handles.labels{ii}).YData-mousepos(3))./Ywidth).^2);
            
            % Show the selected R-peak
            set(handles.temp.R(ii),'XData',handles.graph.R.(handles.labels{ii}).XData(idx),...
                'YData',handles.graph.R.(handles.labels{ii}).YData(idx),...
                'visible','on')
        end
    elseif handles.radiobutton_adjust.Value && length(loop_var) == 1
        for ii = loop_var
            
            % Get the closest R-peak
            [~,idx] = min(((seconds(handles.graph.R.(handles.labels{ii}).XData)-mousepos(1))./Xwidth).^2+...
                ((handles.graph.R.(handles.labels{ii}).YData-mousepos(3))./Ywidth).^2);
            
            % Show the selected R-peak
            set(handles.temp.R(1),'XData',handles.graph.R.(handles.labels{ii}).XData(idx),...
                'YData',handles.graph.R.(handles.labels{ii}).YData(idx),...
                'visible','on')
        end
    else
        % Pre-allocate
        idx = zeros(1,length(loop_var));
        difference = seconds(zeros(1,length(loop_var)));
        for ii = loop_var
            % Get the closest R-peak
            [difference(ii),idx(ii)] = min(((seconds(handles.graph.R.(handles.labels{ii}).XData)-mousepos(1))./Xwidth).^2+...
                ((handles.graph.R.(handles.labels{ii}).YData-mousepos(3))./Ywidth).^2);
        end
        
        % Find the channel
        [~,channel] = min(difference);
        
        % Show the selected R-peak
        set(handles.temp.R(1),'XData',handles.graph.R.(handles.labels{channel}).XData(idx(channel)),...
            'YData',handles.graph.R.(handles.labels{channel}).YData(idx(channel)),...
            'visible','on')
    end
    
else
    % Make the to be deleted R-peak invisible
    set(handles.temp.R,'XData',seconds(NaN))
end

% Update structure
guidata(hObject,handles);


function change_buttonmotionfcn(hObject,eventdata,handles) %#ok
% windowbuttondownfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    % Get the start value
    start = get(handles.temp.R,'XData');
    
    % To prevent trouble in the next step, make start a cell if it is not
    if ~isa(start,'cell')
        start = {start};
    end
    
    % Set the window button motion function
    set(handles.figure_main,'WindowButtonMotionFcn',{@select_R_peaks,handles,start});
end

% Update structure
guidata(hObject,handles);


function select_R_peaks(hObject,eventdata,handles,start) %#ok
% windowbuttonmotionfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    if ~listbox_value % If all leads are displayed
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    % Loop over the selected channels
    for ii = loop_var
        
        % Define start and stop of the search interval
        if seconds(mousepos(1)) < start{ii}
            stop = start{ii};
            start{ii} = seconds(mousepos(1));
        else
            stop = seconds(mousepos(1));
        end
        
        % Find the R-peaks in the search interval
        temp = find(handles.graph.R.(handles.labels{ii}).XData >= start{ii});
        idx = temp(handles.graph.R.(handles.labels{ii}).XData(temp) <= stop);
        
        % Show the selected R-peak(s)
        set(handles.temp.R(ii),'XData',handles.graph.R.(handles.labels{ii}).XData(idx),...
            'YData',handles.graph.R.(handles.labels{ii}).YData(idx),...
            'visible','on')
    end
    
else
    % Make the to be deleted R-peak invisible
    set(handles.temp.R,'XData',seconds(NaN),...
        'YData',handles.data.signal.filtered(1))
end

% Update structure
guidata(hObject,handles);


function delete_R_peak(hObject,eventdata,handles) %#ok
% windowbuttonupfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));
Xwidth = Xlim(2) - Xlim(1);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2) - Ylim(1);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    if ~listbox_value % If all leads are displayed
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    for ii = loop_var
        
        % Get the start
        start = handles.temp.R(ii).XData(1);
        stop = handles.temp.R(ii).XData(end);
        
        % Find the R-peaks on the graph in the search interval
        temp = find(handles.graph.R.(handles.labels{ii}).XData >= start);
        idx = temp(handles.graph.R.(handles.labels{ii}).XData(temp) <= stop);
        
        % Remove this R-peak from the ECG graph
        handles.graph.R.(handles.labels{ii}).XData(idx) = [];
        handles.graph.R.(handles.labels{ii}).YData(idx) = [];
        
        % Remove this RR-interval from the tachogram
        if idx(1) == 1 % First R-peak
            handles.graph.RR.(handles.labels{ii}).XData(idx) = [];
            handles.graph.RR.(handles.labels{ii}).YData(idx) = [];
        elseif idx(end) == length(handles.graph.RR.(handles.labels{ii}).XData) % Last R-peak
            handles.graph.RR.(handles.labels{ii}).XData(min([idx(1) idx(end)-1]):idx(end)-1) = [];
            handles.graph.RR.(handles.labels{ii}).YData(min([idx(1) idx(end)-1]):idx(end)-1) = [];
        else % Every other R-peak
            
            % Adjust the tachogram based on the units selection
            contents = get(handles.popupmenu_units_tachogram,'string');
            switch contents{get(handles.popupmenu_units_tachogram,'Value')}
                case 'HR (bpm)'
                    set(handles.graph.RR.(handles.labels{ii}),'XData',handles.graph.R.(handles.labels{ii}).XData(2:end),...
                        'YData',60./seconds(diff(handles.graph.R.(handles.labels{ii}).XData)));
                otherwise
                    set(handles.graph.RR.(handles.labels{ii}),'XData',handles.graph.R.(handles.labels{ii}).XData(2:end),...
                        'YData',1000*seconds(diff(handles.graph.R.(handles.labels{ii}).XData)));
            end
        end
        
        % Get the new closest R-peak
        [~,idx] = min(((seconds(handles.graph.R.(handles.labels{ii}).XData)-mousepos(1))./Xwidth).^2+...
            ((handles.graph.R.(handles.labels{ii}).YData-mousepos(3))./Ywidth).^2);
        
        % Show the to be deleted R-peak
        set(handles.temp.R(ii),'XData',handles.graph.R.(handles.labels{ii}).XData(idx),...
            'YData',handles.graph.R.(handles.labels{ii}).YData(idx))
    end
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
    
    % Set the window button motion function
    set(handles.figure_main,'WindowButtonMotionFcn',{@select_R_peak,handles},...
        'WindowButtonDownFcn',{@change_buttonmotionfcn,handles});
end

% Update structure
guidata(hObject,handles);

% Add

function drag_new_R_peak(hObject,eventdata,handles) %#ok
% windowbuttonmotionfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));

% Get the y-limits
Ylim = ylim(handles.axes_ECG);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    if mousepos(1) <= Xlim(1) % Past the left limit
        idx = seconds(Xlim(1) + 1/handles.data.fs);
    elseif mousepos(1) >= Xlim(2) % Past the right limit
        idx = seconds(Xlim(2) - 1/handles.data.fs);
    else % Perfect
        idx = seconds(mousepos(1));
    end
    
    % Set up the looping variables
    if ~listbox_value
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    % Show the new R-peak
    for ii = loop_var
        set(handles.temp.R(ii),'XData',idx,...
            'YData',handles.data.signal.filtered(round(seconds(idx)*handles.data.fs),ii),...
            'visible','on')
    end
    
    % Adjust the line
    set(handles.temp.cursor_line,'XData',[idx idx],...
        'YData',Ylim,...
        'visible','on');
else
    % Make the new R-peak(s) invisible
    set(handles.temp.R,'XData',seconds(NaN))
    
    % Make the cursor line(s) invisible
    set(handles.temp.cursor_line,'XData',seconds([NaN NaN]));
end

% Update structure
guidata(hObject,handles);


function define_new_R_peak(hObject,eventdata,handles) %#ok
% windowbuttondownfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));
Xwidth = round((Xlim(2)-Xlim(1))*handles.data.fs);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2)-Ylim(1);

ratio = handles.axes_ECG.Position(3)/handles.axes_ECG.Position(4);

% Check if it is the correct axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    % Get the time variable
    time = seconds((1:size(handles.data.signal.filtered,1))/handles.data.fs);
    
    % Get the x-and y-location
    x_loc = round(mousepos(1)*handles.data.fs); % In samples
    
    % Define the looping variable
    if ~listbox_value
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    for ii = loop_var
        % Create a temporary R-peaks variable in samples
        R_peak_temp = round(seconds(handles.graph.R.(handles.labels{ii}).XData)*handles.data.fs);
        
        % Select a small window (100ms) on each side of the mouse position
        window_width = floor(0.100*handles.data.fs);
        start = max([1 x_loc-window_width]);
        stop = min([x_loc+window_width size(handles.data.signal.filtered,1)]);
        interval = start:stop;
        
        [~,id] = min((((interval-x_loc).*ratio)./Xwidth).^2+...
            ((handles.data.signal.filtered(interval,ii)'-mousepos(3))./Ywidth).^2);
        
        % Get the actual position
        idx = interval(id);
        
        % Locate the new R-peak and add
        id = find(R_peak_temp < idx,1,'last');
        
        if isempty(id) % Before the first R-peak
            R_peak_temp = [idx,R_peak_temp]; %#ok
            
        elseif id == length(R_peak_temp) % After the last R-peak
            R_peak_temp = [R_peak_temp,idx]; %#ok
            
        else
            id = id(end);
            R_peak_temp = [R_peak_temp(1:id),idx,R_peak_temp(id+1:end)];
        end
        
        % Adjust ECG axis
        set(handles.graph.R.(handles.labels{ii}),'XData',time(R_peak_temp),...
            'YData',handles.data.signal.filtered(R_peak_temp,ii),...
            'tag',handles.labels{ii},...
            'UIContextMenu', handles.context_R_peak_options)
        
        % Adjust the tachogram based on the units selection
        contents = get(handles.popupmenu_units_tachogram,'string');
        switch contents{get(handles.popupmenu_units_tachogram,'Value')}
            case 'HR (bpm)'
                set(handles.graph.RR.(handles.labels{ii}),'XData',handles.graph.R.(handles.labels{ii}).XData(2:end),...
                    'YData',60./seconds(diff(handles.graph.R.(handles.labels{ii}).XData)))
            otherwise
                set(handles.graph.RR.(handles.labels{ii}),'XData',handles.graph.R.(handles.labels{ii}).XData(2:end),...
                    'YData',1000*seconds(diff(handles.graph.R.(handles.labels{ii}).XData)))
        end
    end
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end

% Update structure
guidata(hObject,handles);


%% Helper functions adjust R-peak

function define_R_peak(hObject,eventdata,handles) %#ok
% windowbuttondownfcn

% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1;

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = seconds(xlim(handles.axes_ECG));
Xwidth = Xlim(2)-Xlim(1);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2)-Ylim(1);

ratio = handles.axes_ECG.Position(3)/handles.axes_ECG.Position(4);

% Check if there has been a click in the axes
if mousepos(1,1) >= Xlim(1) &&...
        mousepos(1,1) <= Xlim(2) &&...
        mousepos(1,2) >= Ylim(1) &&...
        mousepos(1,2) <= Ylim(2)
    
    if ~listbox_value % If all leads are displayed
        loop_var = 1:handles.data.channels;
    else
        loop_var = listbox_value;
    end
    
    if length(loop_var) == 1
        % Select the closest sample
        [~,handles.temp.to_be_adjusted_sample] = min(((seconds(handles.graph.R.(handles.labels{loop_var}).XData)-mousepos(1))*ratio./Xwidth).^2+...
            ((handles.graph.R.(handles.labels{loop_var}).YData-mousepos(3))./Ywidth).^2);
        
        % Define channel
        channel = loop_var;
    else
        % Pre-allocate
        idx = zeros(1,handles.data.channels);
        difference = seconds(zeros(1,length(loop_var)));
        
        for ii = loop_var
            % Get the closest R-peak
            [difference(ii),idx(ii)] = min(((seconds(handles.graph.R.(handles.labels{ii}).XData)-mousepos(1))*ratio./Xwidth).^2+...
                ((handles.graph.R.(handles.labels{ii}).YData-mousepos(3))./Ywidth).^2);
        end
        
        % Find the channel
        [~,channel] = min(difference);
        
        % Select the closest sample
        handles.temp.to_be_adjusted_sample = idx(channel);
    end
    
    % Set the window button motion function
    set(handles.figure_main,'WindowButtonMotionFcn',{@drag_R_peak,handles,channel});
end

% Update structure
guidata(hObject,handles)


function drag_R_peak(hObject,eventdata,handles,channel) %#ok
% windowbuttonmotionfcn

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = xlim(handles.axes_ECG);
Xwidth = seconds(Xlim(2) - Xlim(1))*handles.data.fs;

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2) - Ylim(1);

ratio = handles.axes_ECG.Position(3)/handles.axes_ECG.Position(4);

% Get the left limit of the search interval
try
    left_lim = handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample-1);
catch
    left_lim = handles.data.start_analysis;
end

% Get the right limit of the search interval
try
    right_lim = handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample+1);
catch
    right_lim = handles.data.start_analysis + handles.data.duration_analysis;
end

if seconds(mousepos(1)) <= left_lim % Past the left limit
    idx = left_lim+seconds(1/handles.data.fs);
elseif seconds(mousepos(1)) >= right_lim % Past the right limit
    idx = right_lim-seconds(1/handles.data.fs);
else % Perfect
    % Get the x-location
    x = round(mousepos(1)*handles.data.fs);
    
    % Adjust the left and right limits
    left_lim = round(seconds(left_lim)*handles.data.fs);
    right_lim = round(seconds(right_lim)*handles.data.fs);
    
    % Select a small window (150ms) on each side of the mouse position
    window_width = floor(0.150*handles.data.fs);
    start = max([1 left_lim x-window_width]);
    stop = min([x+window_width right_lim]);
    interval = start:stop;
    
    [~,id] = min((((interval-x).*ratio)./Xwidth).^2+...
        ((handles.data.signal.filtered(interval,channel)'-mousepos(3))./Ywidth).^2);
    
    idx = seconds(interval(id)/handles.data.fs);
end

% Adjust the temporary R-peak
set(handles.temp.R(1),'XData',idx,...
    'YData',handles.data.signal.filtered(round(seconds(idx)*handles.data.fs),channel),...
    'visible','on')

% Adjust the cursor lines
set(handles.temp.cursor_line,'XData',[idx idx],...
    'YData',Ylim,...
    'visible','on')

% Adjust ECG axis
handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample) = idx;
handles.graph.R.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample) = handles.data.signal.filtered(round(seconds(idx)*handles.data.fs),channel);

% Adjust the tachogram based on the units selection
contents = get(handles.popupmenu_units_tachogram,'string');

if handles.temp.to_be_adjusted_sample == 1
    % YData
    switch contents{get(handles.popupmenu_units_tachogram,'Value')}
        case 'HR (bpm)'
            handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample) = ...
                60./seconds(handles.graph.R.(handles.labels{channel}).XData(2)-handles.graph.R.(handles.labels{channel}).XData(1));
        otherwise
            handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample) = ...
                1000*seconds(handles.graph.R.(handles.labels{channel}).XData(2)-handles.graph.R.(handles.labels{channel}).XData(1));
    end
else
    % XData
    handles.graph.RR.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample-1) = idx;
    
    % YData
    switch contents{get(handles.popupmenu_units_tachogram,'Value')}
        case 'HR (bpm)'
            handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample-1) = ...
                60./seconds(handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample)-handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample-1));
            
            if handles.temp.to_be_adjusted_sample <= length(handles.graph.RR.(handles.labels{channel}).XData)
                handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample) = ...
                    60./seconds(handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample+1)-handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample));
            end
        otherwise
            handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample-1) = ...
                1000*seconds(handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample)-handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample-1));
            
            if handles.temp.to_be_adjusted_sample <= length(handles.graph.RR.(handles.labels{channel}).XData)
                handles.graph.RR.(handles.labels{channel}).YData(handles.temp.to_be_adjusted_sample) = ...
                    1000*seconds(handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample+1)-handles.graph.R.(handles.labels{channel}).XData(handles.temp.to_be_adjusted_sample));
            end
    end
end

% Update structure
guidata(hObject,handles)


function let_R_peak_go(hObject,eventdata,handles)
% windowbuttonup/downfcn

% Remove cursor line
handles.temp.cursor_line.XData = seconds([NaN NaN]);

switch eventdata.EventName
    case 'WindowMouseRelease'
        % Remove the window button motion function
        set(handles.figure_main,'WindowButtonMotionFcn',{@select_R_peak,handles});
        
    case 'WindowMousePress'
        
        % Adjust the model
        handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
end

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Update structure
guidata(hObject,handles);


%% Context axes functions

% --------------------------------------------------------------------
function context_axes_options_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_axes_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current grid situation
set(handles.context_x_grid,'Checked',get(gca,'XGrid'))
set(handles.context_y_grid,'Checked',get(gca,'YGrid'))

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function context_grid_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_x_grid and context_y_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create grids when necessary
switch hObject.Tag
    case 'context_x_grid'
        if strcmp(get(handles.context_x_grid,'Checked'),'on')
            set(handles.context_x_grid,'Checked','off')
            set(gca,'XGrid','off')
        else
            set(handles.context_x_grid,'Checked','on')
            set(gca,'XGrid','on')
        end
    case 'context_y_grid'
        if strcmp(get(handles.context_y_grid,'Checked'),'on')
            set(handles.context_y_grid,'Checked','off')
            set(gca,'YGrid','off')
        else
            set(handles.context_y_grid,'Checked','on')
            set(gca,'YGrid','on')
        end
end

% Update handles structure
guidata(hObject, handles);


%% Context R-peak functions

% --------------------------------------------------------------------
function context_R_peak_options_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_R_peak_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the object
obj = gco;

% Get the tag
tag = get(obj,'Tag');

% Select the closest sample
[~,handles.temp.to_be_adjusted_sample] = min(abs(seconds(mousepos(1))-handles.graph.R.(tag).XData));

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function context_delete_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the object
obj = gco;

% Get the tag
tag = get(obj,'Tag');

% Remove this R-peak from the ECG graph
handles.graph.R.(tag).XData(handles.temp.to_be_adjusted_sample) = [];
handles.graph.R.(tag).YData(handles.temp.to_be_adjusted_sample) = [];

% Remove this RR-interval from the tachogram
if handles.temp.to_be_adjusted_sample == 1 % First R-peak
    handles.graph.RR.(tag).XData(handles.temp.to_be_adjusted_sample) = [];
    handles.graph.RR.(tag).YData(handles.temp.to_be_adjusted_sample) = [];
elseif handles.temp.to_be_adjusted_sample == length(handles.graph.RR.(tag).XData) % Last R-peak
    handles.graph.RR.(tag).XData(min([handles.temp.to_be_adjusted_sample handles.temp.to_be_adjusted_sample-1]):handles.temp.to_be_adjusted_sample-1) = [];
    handles.graph.RR.(tag).YData(min([handles.temp.to_be_adjusted_sample handles.temp.to_be_adjusted_sample-1]):handles.temp.to_be_adjusted_sample-1) = [];
else % Every other R-peak
    
    % Adjust the tachogram based on the units selection
    contents = get(handles.popupmenu_units_tachogram,'string');
    switch contents{get(handles.popupmenu_units_tachogram,'Value')}
        case 'HR (bpm)'
            set(handles.graph.RR.(tag),'XData',handles.graph.R.(tag).XData(2:end),...
                'YData',60./seconds(diff(handles.graph.R.(tag).XData)));
        otherwise
            set(handles.graph.RR.(tag),'XData',handles.graph.R.(tag).XData(2:end),...
                'YData',1000*seconds(diff(handles.graph.R.(tag).XData)));
    end
end

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Adjust the model
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Update structure
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_adjust_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_adjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the object
obj = gco;

% Get the tag
tag = get(obj,'Tag');

% Compare the tag with the labels
listbox_value = find(strcmp(tag,handles.labels));

% Create a temporary R-peak
for ii = 1:handles.data.channels
    handles.temp.R(ii) = plot(handles.axes_ECG,seconds(NaN),handles.data.signal.filtered(1),...
        'oc',...
        'markerfacecolor','c',...
        'linewidth',3);
end

% Set up cursor line
handles.temp.cursor_line = line(handles.axes_ECG,seconds([NaN NaN]), ylim(handles.axes_ECG), ...
    'Color', 'black',...
    'Parent', handles.axes_ECG);

% Set the window button motion function
set(handles.figure_main,'WindowButtonMotionFcn',{@drag_R_peak,handles,listbox_value},...
    'WindowButtonDownFcn',{@let_R_peak_go,handles});

% Update structure
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_converge_to_peak_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_converge_to_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the object
obj = gco;

% Get the tag
tag = get(obj,'Tag');

% Compare the tag with the labels
listbox_value = find(strcmp(tag,handles.labels));

% Define sampling frequency
fs = handles.data.fs;

% Define the signal
signal = handles.data.signal.filtered(:,listbox_value); %#ok

% Define the time
time = seconds((1:size(signal))/fs);
time.Format = 'hh:mm:ss';

% Define a temporary R-peak variable
R_peak_temp = round(fs*seconds(handles.graph.R.(tag).XData(handles.temp.to_be_adjusted_sample)));

% Converge...
if abs(signal(R_peak_temp)) < abs(signal(R_peak_temp-1)) && ...
        abs(signal(R_peak_temp)) < abs(signal(R_peak_temp+1))
    if signal(R_peak_temp) > signal(R_peak_temp-1) && ...
            signal(R_peak_temp) < signal(R_peak_temp+1)
        R_peak_temp = R_peak_temp + 1;
    elseif signal(R_peak_temp) < signal(R_peak_temp-1) && ...
            signal(R_peak_temp) > signal(R_peak_temp+1)
        R_peak_temp = R_peak_temp - 1;
    end
end


% ... forward
while abs(signal(R_peak_temp)) > abs(signal(R_peak_temp-1)) && ...
        abs(signal(R_peak_temp)) < abs(signal(R_peak_temp + 1))
    R_peak_temp = R_peak_temp + 1;
end

% ... backward
while abs(signal(R_peak_temp)) < abs(signal(R_peak_temp-1)) && ...
        abs(signal(R_peak_temp)) > abs(signal(R_peak_temp+1))
    R_peak_temp = R_peak_temp-1;
end

% Store the R-peaks and compute RR-intervals
handles.graph.R.(tag).XData(handles.temp.to_be_adjusted_sample) = time(R_peak_temp);
handles.graph.R.(tag).YData(handles.temp.to_be_adjusted_sample) = signal(R_peak_temp);

% Adjust the model
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Adjust the plots
handles = update_plots(handles);

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Update structure
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_search_extremum_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_search_extremum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the object
obj = gco;

% Get the tag
tag = get(obj,'Tag');

% Compare the tag with the labels
listbox_value = find(strcmp(tag,handles.labels));

% Define sampling frequency
fs = handles.data.fs;

% Define the signal
signal = handles.data.signal.filtered(:,listbox_value); %#ok

% Define the time
time = seconds((1:size(signal))/fs);
time.Format = 'hh:mm:ss';

% Define the window
window = round(0.06*fs);

% Define a temporary R-peak variable
R_peak_temp = round(fs*seconds(handles.graph.R.(tag).XData(handles.temp.to_be_adjusted_sample)));

% Get the samples from the R-peak minus and plus the window
p = signal(max(1,R_peak_temp-window):min(R_peak_temp+window,length(signal)));

% Do the same but for the time locations
t = max(1,R_peak_temp-window):min(R_peak_temp+window,length(signal));

% Do an extremum search based on the selection
switch hObject.Tag
    case 'context_maximum'
        [~,ip] = max(p);
    case 'context_minimum'
        [~,ip] = min(p);
    otherwise
        [~,ip] = max(abs(p));
end

% Adjust the R-peak
R_peak_temp = t(ip(end));

% Store the R-peaks and compute RR-intervals
handles.graph.R.(tag).XData(handles.temp.to_be_adjusted_sample) = time(R_peak_temp);
handles.graph.R.(tag).YData(handles.temp.to_be_adjusted_sample) = signal(R_peak_temp);

% Adjust the model
handles = radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);

% Adjust the plots
handles = update_plots(handles);

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Update structure
guidata(hObject,handles);


%% Picture

% --- Executes on button press in pushbutton_picture.
function pushbutton_picture_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_picture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create a new invisible figure
hfig = figure('Visible','off',...
    'Color','white',...
    'units','norm',...
    'Position',[0 0 1 1]);

% Copy both axes to the new figure
ax(1) = copyobj(handles.axes_ECG,hfig);
ax(2) = copyobj(handles.axes_tachogram,hfig);

% Change the figure and axes units to pixels
set([hfig,ax],'Units','pixels')

% Get the position and tight inset
pos1 = ax(1).Position;
pos2 = ax(2).Position;
ti1 = ax(1).TightInset;
ti2 = ax(2).TightInset;

% Adjust axes horizontal position
ax(1).Position(1) = max([ti1(1) ti2(1)]);
ax(2).Position(1) = max([ti1(1) ti2(1)]);

% Adjust axes vertical position
diff = pos2(2) - ti2(2);
ax(2).Position(2) = ti2(2);
ax(1).Position(2) = pos1(2)-diff;

% Adjust the figure height and width
hfig.Position(3:4) = [pos1(3)+max([ti1(1) ti2(1)])+ti2(3), pos1(2)-pos2(2)+pos1(4)+ti2(2)+ti1(4)];

% Define the name of the picture
[FileName,PathName] = uiputfile({'*.eps','EPS';...
    '*.png','PNG';...
    '*.fig','FIG';...
    '*.jpg','JPG';...
    '*.pdf','PDF'},...
    'Save axes');

% Check if cancel is not pressed
if ~isnumeric(FileName) && ~isnumeric(PathName)
    % Save the picture
    if strcmp(FileName(end-3:end),'.eps')
        set(hfig,'PaperPositionMode','auto')
        print(hfig, fullfile(PathName,FileName), '-depsc ','-painters');
        
    elseif strcmp(FileName(end-3:end),'.pdf')
        set(hfig, 'PaperOrientation', 'landscape');
        print(hfig, fullfile(PathName,FileName), '-dpdf', '-bestfit');
        
    else
        set(hfig,'PaperPositionMode','auto')
        print(hfig, fullfile(PathName,FileName));
    end
end

% Delete figure
clear hfig


%% Link between slider and zoom

% --- Executes on button press in any of the zoom buttons.
function zoom_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to zoom buttons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Source.Tag
    case 'uitoggletool_zoom_in'
        switch hObject.State
            case 'on'
                if ~strcmp(handles.zoom.Enable,'on')
                    % Remove the figure functions
                    set(handles.figure_main,'WindowButtonMotionFcn','',...
                        'WindowButtonDownFcn','',...
                        'WindowButtonUpFcn','');
                end
                
                % Reset both axes
                axes(handles.axes_ECG);
                zoom 'reset'
                axes(handles.axes_tachogram);
                zoom 'reset'
                
                % Set the zoom
                handles.zoom.Direction = 'in';
                handles.zoom.Enable = 'on';
                
            otherwise
                handles.zoom.Enable = 'off';
                
                % Restore the figure functions if necessary
                if handles.radiobutton_add.Value+handles.radiobutton_adjust.Value+handles.radiobutton_delete.Value == 1
                    radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
                elseif handles.togglebutton_analysis_period.Value
                    if ~strcmp(handles.temp.start.Visible,'on')
                        set(handles.figure_main,'WindowButtonMotionFcn',{@drag_start_line,handles},...
                            'WindowButtonDownFcn',{@define_start_line,handles},...
                            'WindowButtonUpFcn',{@define_stop_line,handles});
                    else
                        set(handles.figure_main,'WindowButtonDownFcn',{@select_line,handles},...
                            'WindowButtonUpFcn',{@define_stop_line,handles});
                    end
                end
        end
    case 'uitoggletool_zoom_out'
        switch hObject.State
            case 'on'
                if ~strcmp(handles.zoom.Enable,'on')
                    % Remove the figure functions
                    set(handles.figure_main,'WindowButtonMotionFcn','',...
                        'WindowButtonDownFcn','',...
                        'WindowButtonUpFcn','');
                end
                
                % Set the zoom
                handles.zoom.Direction = 'out';
                handles.zoom.Enable = 'on';
                
            otherwise
                handles.zoom.Enable = 'off';
                
                % Restore the figure functions if necessary
                if handles.radiobutton_add.Value+handles.radiobutton_adjust.Value+handles.radiobutton_delete.Value == 1
                    radiobutton_add_adjust_delete_Callback(hObject, eventdata, handles);
                elseif handles.togglebutton_analysis_period.Value
                    if ~strcmp(handles.temp.start.Visible,'on')
                        set(handles.figure_main,'WindowButtonMotionFcn',{@drag_start_line,handles},...
                            'WindowButtonDownFcn',{@define_start_line,handles},...
                            'WindowButtonUpFcn',{@define_stop_line,handles});
                    else
                        set(handles.figure_main,'WindowButtonDownFcn',{@select_line,handles},...
                            'WindowButtonUpFcn',{@define_stop_line,handles});
                    end
                end
        end
end

% Update handles structure
guidata(hObject, handles);


function zoom_postcallback(hObject, eventdata)

% Get the handles
handles = guidata(eventdata.Axes);

% Get the start of the analysis period in duration values
an_start = handles.data.start_analysis;

% Get the width of the analysis window in duration values
an_width = handles.data.duration_analysis;

% Get the new limits
newLim = eventdata.Axes.XLim;

% Get the new range
range = newLim(2)-newLim(1);

% Check the new range
if range < seconds(1)
    
    % Set the range
    range = seconds(1);
    range.Format = 'hh:mm:ss';
    
    if newLim(1) + range > an_start + an_width
        newLim(1) = an_start + an_width - seconds(1);
        newLim(2) = an_start + an_width;
    else
        newLim(2) = newLim(1) + range;
    end
    
    % Change the limits
    set(eventdata.Axes,'XLim',newLim)
end

% Adjust the range
set(handles.edit_range,'string',char(range));

% Adjust the data
handles.data.range = range;

if ~isempty(handles.data.signal.filtered)
    
    % Get a new maximum for the slider
    new_max = round(handles.data.fs*seconds(an_width-range));
    
    % Get the new step size
    n = round(handles.data.fs*seconds(range))/new_max;
    
    % Adjust the slider
    set(handles.slider_range,'max',new_max,...
        'sliderstep',[min([1 n]) min([1 2*n])],...
        'string','zoom',...
        'value',round(handles.data.fs*seconds(newLim(1)-an_start)))
    
    % Adjust the string
    set(handles.slider_range,'string','normal')
    
    % Adjust the slider's enability
    if ~strcmp(char(handles.data.duration_analysis),char(range))
        set(handles.slider_range,'Enable','on')
    else
        set(handles.slider_range,'Enable','off')
    end
    
    % Adjust y-limits of the tachogram if the ECG axes is clicked and if
    % the y-values are not fixed
    if strcmp(eventdata.Axes.Tag,'axes_ECG')
        
        handles = adjust_y_limits(handles,{'axes_tachogram'});
        
    elseif strcmp(eventdata.Axes.Tag,'axes_tachogram') % Adjust y-limits of the ECG graph
        
        handles = adjust_y_limits(handles,{'axes_ECG'});
    end
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function context_lead_options_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_lead_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function context_add_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the mouse position
mousepos = get(handles.axes_ECG,'CurrentPoint');

% Get the x-limits
Xlim = xlim(handles.axes_ECG);
Xwidth = round(seconds(Xlim(2) - Xlim(1))*handles.data.fs);

% Get the y-limits
Ylim = ylim(handles.axes_ECG);
Ywidth = Ylim(2) - Ylim(1);

% Get the ratio
ratio = handles.axes_ECG.Position(3)/handles.axes_ECG.Position(4);

% Get the selected object
current_obj = gco;

% Get the tag
tag = current_obj.Tag;

% Get the selected lead
switch tag
    case handles.labels(1)
        lead = 1;
    case handles.labels(2)
        lead = 2;
    case handles.labels(3)
        lead = 3;
    case handles.labels(4)
        lead = 4;
    case handles.labels(5)
        lead = 5;
    case handles.labels(6)
        lead = 6;
    case handles.labels(7)
        lead = 7;
    case handles.labels(8)
        lead = 8;
    case handles.labels(9)
        lead = 9;
    case handles.labels(10)
        lead = 10;
    case handles.labels(11)
        lead = 11;
    case handles.labels(12)
        lead = 12;
end

% Get the time variable
time = seconds((1:size(handles.data.signal.filtered,1))/handles.data.fs);

% Get the x-location in samples
x = round(mousepos(1)*handles.data.fs);

% Select a small window (150ms) on each side of the mouse position
window_width = floor(0.150*handles.data.fs);
start = max([1 x-window_width]);
stop = min([x+window_width size(handles.data.signal.filtered,1)]);
interval = start:stop;

[~,id] = min((((interval-x).*ratio)./Xwidth).^2+...
    ((handles.data.signal.filtered(interval,lead)'-mousepos(3))./Ywidth).^2);

idx = interval(id);

% Create a temporary R-peaks variable in samples
R_peak_temp = round(seconds(handles.graph.R.(handles.labels{lead}).XData)*handles.data.fs);

% Locate the new R-peak and add
id = find(R_peak_temp < idx,1,'last');

if isempty(id) % Before the first R-peak
    R_peak_temp = [idx,R_peak_temp];
    
elseif id == length(R_peak_temp) % After the last R-peak
    R_peak_temp = [R_peak_temp,idx];
else
    id = id(end);
    R_peak_temp = [R_peak_temp(1:id),idx,R_peak_temp(id+1:end)];
end

% Adjust ECG axis
set(handles.graph.R.(handles.labels{lead}),'XData',time(R_peak_temp),...
    'YData',handles.data.signal.filtered(R_peak_temp,lead),...
    'tag',handles.labels{lead},...
    'UIContextMenu', handles.context_R_peak_options)

% Adjust the tachogram based on the units selection
contents = get(handles.popupmenu_units_tachogram,'string');
switch contents{get(handles.popupmenu_units_tachogram,'Value')}
    case 'HR (bpm)'
        set(handles.graph.RR.(handles.labels{lead}),'XData',handles.graph.R.(handles.labels{lead}).XData(2:end),...
            'YData',60./seconds(diff(handles.graph.R.(handles.labels{lead}).XData)))
    otherwise
        set(handles.graph.RR.(handles.labels{lead}),'XData',handles.graph.R.(handles.labels{lead}).XData(2:end),...
            'YData',1000*seconds(diff(handles.graph.R.(handles.labels{lead}).XData)))
end

% Adjust the y-limits
handles = adjust_y_limits(handles);

% Update structure
guidata(hObject,handles);


% --------------------------------------------------------------------
function context_change_color_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_change_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the selected object
current_obj = gco;

% Get the tag
tag = current_obj.Tag;

% Get the color
current_color = get(current_obj,'Color');

% Select a new color
new_color = uisetcolor(current_color,'Select a new color');

% Get the selected lead
switch tag
    case handles.labels(1)
        idx = 1;
    case handles.labels(2)
        idx = 2;
    case handles.labels(3)
        idx = 3;
end

% Adjust the color variable
handles.color(idx,:) = new_color;

% Adjust the listbox

old_string = handles.listbox_channels.String{idx+1};
name = old_string(find(old_string=='"',1,'last')+2:end-14);
new_colorStr = sprintf('%d,',int16(255*new_color));
handles.listbox_channels.String{idx+1} = ['<HTML><FONT color="rgb(' new_colorStr ')">' name '</Font></html>'];

% Update the plots
handles = update_plots(handles);

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function context_delete_lead_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to context_delete_lead (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the selected object
current_obj = gco;

% Get the tag
tag = current_obj.Tag;

% Delete selected lead plot and the original lead plot
delete(handles.graph.signal.filtered.(tag))
delete(handles.graph.signal.original.(tag))

% Delete R-peaks if present
if ~isempty(handles.graph.R)
    delete(handles.graph.R.(tag))
    delete(handles.graph.RR.(tag))
end

% Get the to be deleted lead
switch tag
    case handles.labels(1)
        idx=1;
    case handles.labels(2)
        idx=2;
    case handles.labels(3)
        idx=3;
end

% Delete the lead
handles.data.signal.filtered(:,idx) = [];
handles.data.signal.original(:,idx) = [];

% Adjust the lead's label
handles.labels(idx) = [];

% Subtract one value of the nr of channels
handles.data.channels = handles.data.channels-1;

% Delete the color
handles.color(idx,:) = [];

% Adjust the listbox
handles.listbox_channels.String(idx+1) = [];

set(handles.listbox_channels,'value',1,...
    'enable','on')

% Check if there is still a signal present. If not, disable everything.
if ~handles.data.channels
    set(findall(handles.figure_main, '-property', 'enable'), 'enable', 'off')
    
    % Enable the menu tools
    set(findall(handles.figure_main,'Type','uimenu'),'enable','on')
    
else
    % Update the plots
    handles = update_plots(handles);
    
    % Adjust the y-limits
    handles = adjust_y_limits(handles);
end


% Update handles structure
guidata(hObject, handles);


% --- Executes when figure_main is resized.
function figure_main_SizeChangedFcn(hObject, eventdata, handles) %#ok
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if ~handles.uiLimitSet
        LimitFigSize(handles.figure_main, 'min', [1100, 600])
    end
catch
end


% --- Executes when user attempts to close figure_main.
function figure_main_CloseRequestFcn(hObject, eventdata, handles) %#ok
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% Check if a signal has been loaded and if the results are already saved
choice = 'Yes';
if ~handles.check.save && ~isempty(handles.data.signal.original)
    choice = questdlg('Current analysis has not been saved. Are you sure you want to quit?', ...
        'Quit', ...
        'Yes','No','Cancel','Yes');
end

% Next step depends on the answer of the previous question
switch choice
    case 'Yes'
        % Close the preferences figure
        try
            if isvalid(handles.preferences.Fig)
                % Delete the app
                close(handles.preferences.Fig)
            end
        catch
        end
        
        delete(handles.figure_main)
end
