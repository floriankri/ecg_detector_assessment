function varargout = Preferences_GUIDE(varargin)
% PREFERENCES_GUIDE MATLAB code for Preferences_GUIDE.fig
%      PREFERENCES_GUIDE, by itself, creates a new PREFERENCES_GUIDE or raises the existing
%      singleton*.
%
%      H = PREFERENCES_GUIDE returns the handle to a new PREFERENCES_GUIDE or the handle to
%      the existing singleton*.
%
%      PREFERENCES_GUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREFERENCES_GUIDE.M with the given input arguments.
%
%      PREFERENCES_GUIDE('Property','Value',...) creates a new PREFERENCES_GUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Preferences_GUIDE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Preferences_GUIDE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Preferences_GUIDE

% Last Modified by GUIDE v2.5 24-May-2019 13:22:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Preferences_GUIDE_OpeningFcn, ...
    'gui_OutputFcn',  @Preferences_GUIDE_OutputFcn, ...
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

% --- Executes just before Preferences_GUIDE is made visible.
function Preferences_GUIDE_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Preferences_GUIDE (see VARARGIN)

% Set the background color 
handles.color = [0.07 0.72 1];

% Store the parameters
parameters = varargin{1,1};

% Set up the Welch tab
if isempty(parameters.Welch{1})
    % Disable the settings
    handles.edit_window_width.Enable = 'off';
    handles.text_seconds.Enable = 'off';
    
    % Select the default button
    handles.radiobutton_default_window_width.Value = 1;
else
    % Enable the settings
    handles.edit_window_width.Enable = 'on';
    handles.text_seconds.Enable = 'on';
    
    % Select the other button
    handles.radiobutton_window_width.Value = 1;
    
    % Set the window width
    handles.edit_window_width.String = num2str(parameters.Welch{1});
end

if isempty(parameters.Welch{2})
    % Disable the settings
    handles.edit_overlap.Enable = 'off';
    handles.text_percentage.Enable = 'off';
    
    % Select the default button
    handles.radiobutton_default_overlap.Value = 1;
else
    % Enable the settings
    handles.edit_overlap.Enable = 'on';
    handles.text_percentage.Enable = 'on';
    
    % Select the other button
    handles.radiobutton_overlap.Value = 1;
    
    % Set the overlap
    handles.edit_overlap.Value = num2str(parameters.Welch{2});
end

if isempty(parameters.Welch{3})
    % Disable the settings
    handles.edit_nfft.Enable = 'off';
    handles.text_seconds2.Enable = 'off';
    
    % Select the default button
    handles.radiobutton_default_nfft.Value = 1;
else
    % Enable the settings
    handles.edit_nfft.Enable = 'on';
    handles.text_seconds2.Enable = 'on';
    
    % Select the other button
    handles.radiobutton_nfft.Value = 1;
    
    % Set the nfft
    handles.edit_nfft.Value = num2str(parameters.Welch{3});
end

% Set up the Filter tab
handles.edit_high_pass_order.String = num2str(parameters.Filter{1});
handles.edit_low_pass_order.String = num2str(parameters.Filter{2});

% Set up the R-peak detection tab
handles.edit_envelope_size.String = num2str(parameters.Rpeak{1});
handles.edit_average_hr.String = num2str(parameters.Rpeak{2});
handles.radiobutton_post_yes.Value = parameters.Rpeak{3};
handles.radiobutton_post_no.Value = ~parameters.Rpeak{3};
handles.radiobutton_ectopic_yes.Value = parameters.Rpeak{4};
handles.radiobutton_ectopic_no.Value = ~parameters.Rpeak{4};

% Store the parameters
handles.Parameters = parameters;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Preferences_GUIDE wait for user response (see UIRESUME)
uiwait(handles.figure_main);


% --- Outputs from this function are returned to the command line.
function varargout = Preferences_GUIDE_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set default output
varargout{1} = [];

% If figure is not closed or canceled, set varargout equal to the
% parameters
if ~isempty(handles)
    varargout{1} = handles.Parameters;
    
    % Delete figure
    delete(handles.figure_main)
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_pre_processing.
function text_pre_processing_ButtonDownFcn(hObject, eventdata, handles) %#ok
% hObject    handle to text_pre_processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adjust the background
set(findall(handles.uipanel_selection,'Style','text'),'BackgroundColor',[1 1 1])
handles.text_pre_processing.BackgroundColor = handles.color;

% Show the panel
handles.uipanel_empty.Visible = 'on';

% Make the other panels invisible
handles.uipanel_welch_periodogram.Visible = 'off';
handles.uipanel_high_pass.Visible = 'off';
handles.uipanel_low_pass.Visible = 'off';
handles.uipanel_peak_detection.Visible = 'off';

% Disable the default button
handles.pushbutton_default.Enable = 'off';

switch handles.text_pre_processing.Value
    case 0
        % Adjust the text
        handles.text_pre_processing.String = '<  Pre-processing';
        
        % Show the subtabs
        handles.text_welch_periodogram.Visible = 'on';
        handles.text_filter.Visible = 'on';
        
        % Move the R-peaks tab
        handles.text_r_peak_detection.Position(2) = 13.25;
        
        % Adjust the value
        handles.text_pre_processing.Value = 1;
    case 1
        % Adjust the text
        handles.text_pre_processing.String = '>  Pre-processing';
        
        % Show the subtabs
        handles.text_welch_periodogram.Visible = 'off';
        handles.text_filter.Visible = 'off';
        
        % Move the R-peaks tab
        handles.text_r_peak_detection.Position(2) = 16;
        
        % Adjust the value
        handles.text_pre_processing.Value = 0;
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_welch_periodogram.
function text_welch_periodogram_ButtonDownFcn(hObject, eventdata, handles) %#ok
% hObject    handle to text_welch_periodogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adjust the background
set(findall(handles.uipanel_selection,'Style','text'),'BackgroundColor',[1 1 1])
handles.text_welch_periodogram.BackgroundColor = handles.color;

% Show the panel
handles.uipanel_welch_periodogram.Visible = 'on';

% Make the other panels invisible
handles.uipanel_empty.Visible = 'off';
handles.uipanel_high_pass.Visible = 'off';
handles.uipanel_low_pass.Visible = 'off';
handles.uipanel_peak_detection.Visible = 'off';

% Enable the default button
handles.pushbutton_default.Enable = 'on';


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_filter.
function text_filter_ButtonDownFcn(hObject, eventdata, handles) %#ok
% hObject    handle to text_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adjust the background
set(findall(handles.uipanel_selection,'Style','text'),'BackgroundColor',[1 1 1])
handles.text_filter.BackgroundColor = handles.color;

% Show the panels
handles.uipanel_high_pass.Visible = 'on';
handles.uipanel_low_pass.Visible = 'on';

% Make the other panels invisible
handles.uipanel_empty.Visible = 'off';
handles.uipanel_welch_periodogram.Visible = 'off';
handles.uipanel_peak_detection.Visible = 'off';

% Enable the default button
handles.pushbutton_default.Enable = 'on';


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_r_peak_detection.
function text_r_peak_detection_ButtonDownFcn(hObject, eventdata, handles) %#ok
% hObject    handle to text_r_peak_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adjust the background
set(findall(handles.uipanel_selection,'Style','text'),'BackgroundColor',[1 1 1])
handles.text_r_peak_detection.BackgroundColor = handles.color;

% Show the panels
handles.uipanel_peak_detection.Visible = 'on';

% Make the other panels invisible
handles.uipanel_empty.Visible = 'off';
handles.uipanel_welch_periodogram.Visible = 'off';
handles.uipanel_high_pass.Visible = 'off';
handles.uipanel_low_pass.Visible = 'off';

% Enable the default button
handles.pushbutton_default.Enable = 'on';


% --- Executes when selected object is changed in uibuttongroup_window_width.
function uibuttongroup_window_width_SelectionChangedFcn(hObject, eventdata, handles) %#ok
% hObject    handle to the selected object in uibuttongroup_window_width
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.radiobutton_default_window_width.Value
    case 1
        handles.edit_window_width.Enable = 'off';
        handles.text_seconds.Enable = 'off';
    otherwise
        handles.edit_window_width.Enable = 'on';
        handles.text_seconds.Enable = 'on';
end


function edit_window_width_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_window_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_window_width as text
%        str2double(get(hObject,'String')) returns contents of edit_window_width as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The window width should be a number.')
    
    % Reset to default
    handles.radiobutton_default_window_width.Value = 1;
    if ~isempty(handles.Parameters.Welch{1})
        handles.edit_window_width.String = num2str(handles.Parameters.Welch{1});
    else
        handles.edit_window_width.String = '1';
    end
elseif value <= 0
    % Show a warning dialog
    warndlg('The window width should be larger than zero.')
    
    % Reset to default
    handles.radiobutton_default_window_width.Value = 1;
    if ~isempty(handles.Parameters.Welch{1})
        handles.edit_window_width.String = num2str(handles.Parameters.Welch{1});
    else
        handles.edit_window_width.String = '1';
    end
end

% Run the selection changed function
uibuttongroup_window_width_SelectionChangedFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_window_width_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_window_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup_overlap.
function uibuttongroup_overlap_SelectionChangedFcn(hObject, eventdata, handles) %#ok
% hObject    handle to the selected object in uibuttongroup_overlap
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.radiobutton_default_overlap.Value
    case 1
        handles.edit_overlap.Enable = 'off';
        handles.text_percentage.Enable = 'off';
    otherwise
        handles.edit_overlap.Enable = 'on';
        handles.text_percentage.Enable = 'on';
end


function edit_overlap_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_overlap as text
%        str2double(get(hObject,'String')) returns contents of edit_overlap as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The overlap should be a number.')
    
    % Reset to default
    handles.radiobutton_default_overlap.Value = 1;
    if ~isempty(handles.Parameters.Welch{2})
        handles.edit_overlap.String = num2str(handles.Parameters.Welch{2});
    else
        handles.edit_overlap.String = '50';
    end
elseif value < 1 || value > 100
    % Show a warning dialog
    warndlg('The overlap should be between 1 and 100%')
    
    % Reset to default
    handles.radiobutton_default_overlap.Value = 1;
    if ~isempty(handles.Parameters.Welch{2})
        handles.edit_overlap.String = num2str(handles.Parameters.Welch{2});
    else
        handles.edit_overlap.String = '50';
    end
end

% Run the selection changed function
uibuttongroup_overlap_SelectionChangedFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_overlap_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup_nfft.
function uibuttongroup_nfft_SelectionChangedFcn(hObject, eventdata, handles)%#ok
% hObject    handle to the selected object in uibuttongroup_nfft
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.radiobutton_default_nfft.Value
    case 1
        handles.edit_nfft.Enable = 'off';
        handles.text_seconds2.Enable = 'off';
    otherwise
        handles.edit_nfft.Enable = 'on';
        handles.text_seconds2.Enable = 'on';
end


function edit_nfft_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nfft as text
%        str2double(get(hObject,'String')) returns contents of edit_nfft as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The nfft should be a number.')
    
    % Reset to default
    handles.radiobutton_default_nfft.Value = 1;
    if ~isempty(handles.Parameters.Welch{3})
        handles.edit_nfft.String = num2str(handles.Parameters.Welch{3});
    else
        handles.edit_nfft.String = '1';
    end
elseif value <= 0
    % Show a warning dialog
    warndlg('The nfft should be larger than zero.')
    
    % Reset to default
    handles.radiobutton_default_overlap.Value = 1;
    if ~isempty(handles.Parameters.Welch{3})
        handles.edit_nfft.String = num2str(handles.Parameters.Welch{3});
    else
        handles.edit_nfft.String = '1';
    end
end

% Run the selection changed function
uibuttongroup_nfft_SelectionChangedFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_nfft_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_high_pass_order_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_high_pass_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_high_pass_order as text
%        str2double(get(hObject,'String')) returns contents of edit_high_pass_order as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The order should be a number.')
    
    % Reset to default
    handles.edit_high_pass_order.String = num2str(handles.Parameters.Filter{1});
elseif value <= 0
    % Show a warning dialog
    warndlg('The order should be larger than zero.')
    
    % Reset to default
    handles.edit_high_pass_order.String = num2str(handles.Parameters.Filter{1});
end


% --- Executes during object creation, after setting all properties.
function edit_high_pass_order_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_high_pass_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_low_pass_order_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_low_pass_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_low_pass_order as text
%        str2double(get(hObject,'String')) returns contents of edit_low_pass_order as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The order should be a number.')
    
    % Reset to default
    handles.edit_low_pass_order.String = num2str(handles.Parameters.Filter{2});
elseif value <= 0
    % Show a warning dialog
    warndlg('The order should be larger than zero.')
    
    % Reset to default
    handles.edit_low_pass_order.String = num2str(handles.Parameters.Filter{2});
end


% --- Executes during object creation, after setting all properties.
function edit_low_pass_order_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_low_pass_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_envelope_size_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_envelope_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_envelope_size as text
%        str2double(get(hObject,'String')) returns contents of edit_envelope_size as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The envelope size should be a number.')
    
    % Reset to default
    handles.edit_envelope_size.String = num2str(handles.Parameters.Rpeak{1});
elseif value < 1 || value > 600
    % Show a warning dialog
    warndlg('The envelope size should be between 1 and 600 ms.')
    
    % Reset to default
    handles.edit_envelope_size.String = num2str(handles.Parameters.Rpeak{1});
end

% --- Executes during object creation, after setting all properties.
function edit_envelope_size_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_envelope_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_average_hr_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_average_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_average_hr as text
%        str2double(get(hObject,'String')) returns contents of edit_average_hr as a double

% Get the value
value = str2double(get(hObject,'String'));

if isnan(value)
    % Show a warning dialog
    warndlg('The envelope size should be a number.')
    
    % Reset to default
    handles.edit_average_hr.String = num2str(handles.Parameters.Rpeak{2});
elseif value < 1 || value > 600
    % Show a warning dialog
    warndlg('The average HR should be between 30 and 200 ms.')
    
    % Reset to default
    handles.edit_average_hr.String = num2str(handles.Parameters.Rpeak{2});
end


% --- Executes during object creation, after setting all properties.
function edit_average_hr_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_average_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Store the Welch parameters
if handles.radiobutton_default_window_width.Value == 0
    handles.Parameters.Welch{1} = str2double(handles.edit_window_width.String);
else
    handles.Parameters.Welch{1} = [];
end
if handles.radiobutton_default_overlap.Value == 0
    handles.Parameters.Welch{2} = str2double(handles.edit_overlap.String);
else
    handles.Parameters.Welch{2} = [];
end
if handles.radiobutton_default_nfft.Value == 0
    handles.Parameters.Welch{3} = str2double(handles.edit_nfft.String);
else
    handles.Parameters.Welch{3} = [];
end

% Store the Filter parameters
handles.Parameters.Filter{1} = str2double(handles.edit_high_pass_order.String);
handles.Parameters.Filter{2} = str2double(handles.edit_low_pass_order.String);

% Store the R-peak parameters
handles.Parameters.Rpeak{1} = str2double(handles.edit_envelope_size.String);
handles.Parameters.Rpeak{2} = str2double(handles.edit_average_hr.String);
handles.Parameters.Rpeak{3} = handles.radiobutton_post_yes.Value;
handles.Parameters.Rpeak{4} = handles.radiobutton_ectopic_yes.Value;

% Update handles structure
guidata(hObject, handles);

% Resume GUI
uiresume(handles.figure_main)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Resume GUI
uiresume(handles.figure_main)


% --- Executes on button press in pushbutton_default.
function pushbutton_default_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the selected panel
panel = findall(findall(handles.uipanel_selection,'Style','text'),'BackgroundColor',handles.color);

switch panel.Tag
    case 'text_welch_periodogram'
        % Disable the settings
        handles.edit_window_width.Enable = 'off';
        handles.text_seconds.Enable = 'off';
        handles.edit_overlap.Enable = 'off';
        handles.text_percentage.Enable = 'off';
        handles.edit_nfft.Enable = 'off';
        handles.text_seconds2.Enable = 'off';
        
        % Select the default buttons
        handles.radiobutton_default_window_width.Value = 1;
        handles.radiobutton_default_overlap.Value = 1;
        handles.radiobutton_default_nfft.Value = 1;
    case 'text_filter'
        handles.edit_high_pass_order.String = '2';
        handles.edit_low_pass_order.String = '4';
    case 'text_r_peak_detection'
        handles.edit_envelope_size.String = '300';
        handles.edit_average_hr.String = '100';
        handles.radiobutton_post_yes.Value = 1;
        handles.radiobutton_ectopic_yes.Value = 0;
end
