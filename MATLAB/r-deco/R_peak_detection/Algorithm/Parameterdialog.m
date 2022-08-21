function varargout = Parameterdialog(varargin)
% PARAMETERDIALOG MATLAB code for Parameterdialog.fig
%      PARAMETERDIALOG, by itself, creates a new PARAMETERDIALOG or raises the existing
%      singleton*.
%
%      H = PARAMETERDIALOG returns the handle to a new PARAMETERDIALOG or the handle to
%      the existing singleton*.
%
%      PARAMETERDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERDIALOG.M with the given input arguments.
%
%      PARAMETERDIALOG('Property','Value',...) creates a new PARAMETERDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Parameterdialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Parameterdialog_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help Parameterdialog

% Last Modified by GUIDE v2.5 23-Nov-2018 14:42:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Parameterdialog_OpeningFcn, ...
                   'gui_OutputFcn',  @Parameterdialog_OutputFcn, ...
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


% --- Executes just before Parameterdialog is made visible.
function Parameterdialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Parameterdialog (see VARARGIN)

% Assign the parameters
handles.parameters = varargin{1};

% Get the sampling frequency
handles.fs = varargin{3};

% Get ten seconds of the first lead signal
temp_signal = varargin{2};

if length(temp_signal(:,1))/handles.fs > 10
    % Get the signal
    handles.signal = temp_signal(1:handles.fs*10,1);

    % Get the time
    handles.time = seconds((1:10*handles.fs)/handles.fs);
else
    % Get the signal
    handles.signal = temp_signal(:,1);
    
    % Get the time
    handles.time = seconds((1:length(handles.signal))/handles.fs);
end

% Plot something random
plot(handles.axes_ECG,seconds([0 60]),[0 1],...
    'visible','off');

% Adjust the xtickformat
xtickformat(handles.axes_ECG,'hh:mm:ss')

% Define the axis labels
xlabel(handles.axes_ECG,{'Time (hh:mm:ss)'})

% Remove y-axes
set(handles.axes_ECG,'YTick', [])

% Update the GUI
update_gui(hObject, eventdata, handles)

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Parameterdialog wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Parameterdialog_OutputFcn(hObject, eventdata, handles) %#ok 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set default output
varargout{1} = [];

% If figure is not closed or canceled, set varargout equal to the
% parameters
if ~isempty(handles)
    varargout{1} = handles.parameters;
    
    % Delete figure
    delete(handles.figure1)
end


function update_gui(hObject, eventdata, handles)
% hObject    handle to selected uicontrol
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set up the envelope
set(handles.edit_envelope_size,'string',handles.parameters{1})
pushbutton_apply_Callback(hObject, eventdata,handles)

% Set up the average heart rate
set(handles.edit_average_heart_rate,'string',handles.parameters{2})

% Set up the post processing
if handles.parameters{3}
    set(handles.radiobutton_yes_post,'Value',1)
    set(handles.radiobutton_no_post,'Value',0)
else
    set(handles.radiobutton_yes_post,'Value',0)
    set(handles.radiobutton_no_post,'Value',1)
end

% Set up ectopics
if handles.parameters{4}
    set(handles.radiobutton_yes_ectopic,'Value',1)
    set(handles.radiobutton_no_ectopic,'Value',0)
else
    set(handles.radiobutton_yes_ectopic,'Value',0)
    set(handles.radiobutton_no_ectopic,'Value',1)
end

% Set up inverted
if handles.parameters{5}
    set(handles.radiobutton_yes_inverted,'Value',1)
    set(handles.radiobutton_no_inverted,'Value',0)
else
    set(handles.radiobutton_yes_inverted,'Value',0)
    set(handles.radiobutton_no_inverted,'Value',1)
end

% Update handles structure
guidata(hObject, handles);


function edit_envelope_size_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_envelope_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_envelope_size as text
%        str2double(get(hObject,'String')) returns contents of edit_envelope_size as a double


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


% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temp = str2double(get(handles.edit_envelope_size,'String'));

if isnan(temp) 
    errordlg('The envelope size should be numerical.','Format Error')

    set(handles.edit_envelope_size,'string',handles.parameters{1})
elseif temp <= 0
    errordlg('The envelope size should be positive.','Format Error')

    set(handles.edit_envelope_size,'string',handles.parameters{1})
elseif temp <= 1
    warndlg('The envelope size should be larger than 1 ms.','Warning')

    set(handles.edit_envelope_size,'string',handles.parameters{1})
elseif temp > 600 && temp ~= handles.parameters{1}
    choice = questdlg('Envelope sizes larger than 600 ms could lead to bad R-peak detection. Are you sure you want to continue?',...
        'Warning',...
        'Yes','No','Cancel','No');

    switch choice
        case 'Yes'
            handles.parameters{1} = temp;
        case {'No', 'Cancel', '' }
            set(handles.edit_envelope_size,'string',handles.parameters{1})
    end
else
    handles.parameters{1} = temp;
end

% Set envelope in samples
env = round(handles.fs*handles.parameters{1}/1000);

% Upper envelope
up_env = env_secant(seconds(handles.time),handles.signal,env,'top');

% Lower envelope
low_env = env_secant(seconds(handles.time),handles.signal,env,'bottom');

% Enveloped segment
env_signal = up_env-low_env;

% Clear axes
cla(handles.axes_ECG)

% Plot the signal
plot(handles.axes_ECG,handles.time,handles.signal)

% Plot the enveloped signal
plot(handles.axes_ECG,handles.time,env_signal)

% Set the limits
xlim(handles.axes_ECG,seconds([0 10]))

try
    ylim(handles.axes_ECG,[(min([env_signal'; handles.signal])-std([env_signal'; handles.signal])) (max([env_signal'; handles.signal])+2*std([env_signal'; handles.signal]))])
catch
end

% Set the legend
legend(handles.axes_ECG,'ECG','Envelope')

% Update handles structure
guidata(hObject, handles);


function edit_average_heart_rate_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_average_heart_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_average_heart_rate as text
%        str2double(get(hObject,'String')) returns contents of edit_average_heart_rate as a double

temp = str2double(get(handles.edit_average_heart_rate,'String'));

if isnan(temp) || temp <= 20 || temp >= 250

    set(handles.edit_average_heart_rate,'string',handles.parameters{2})
else
    handles.parameters{2} = temp;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_average_heart_rate_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_average_heart_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_yes_post.
function radiobutton_yes_post_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_yes_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_yes_post
if get(hObject,'Value') == 1
    set(handles.radiobutton_no_post,'Value',0)
    handles.parameters{3} = 1;
else
    set(handles.radiobutton_no_post,'Value',1)
    handles.parameters{3} = 0;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in radiobutton_no_post.
function radiobutton_no_post_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_no_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_no_post
if get(hObject,'Value') == 1
    set(handles.radiobutton_yes_post,'Value',0)
    handles.parameters{3} = 0;
else
    set(handles.radiobutton_yes_post,'Value',1)
    handles.parameters{3} = 1;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in radiobutton_yes_ectopic and radiobutton_no_ectopic.
function radiobutton_ectopic_yes_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_yes_ectopic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_yes_ectopic
if get(hObject,'Value') == 1
    set(handles.radiobutton_no_ectopic,'Value',0)
    handles.parameters{4} = 1;
else
    set(handles.radiobutton_no_ectopic,'Value',1)
    handles.parameters{4} = 0;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in radiobutton_no_ectopic.
function radiobutton_ectopic_no_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_yes_ectopic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_no_ectopic
if get(hObject,'Value') == 1
    set(handles.radiobutton_yes_ectopic,'Value',0)
    handles.parameters{4} = 0;
else
    set(handles.radiobutton_yes_ectopic,'Value',1)
    handles.parameters{4} = 1;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in radiobutton_yes_inverted.
function radiobutton_yes_inverted_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_yes_inverted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_yes_inverted
if get(hObject,'Value') == 1
    set(handles.radiobutton_no_inverted,'Value',0)
    handles.parameters{5} = 1;
else
    set(handles.radiobutton_no_inverted,'Value',1)
    handles.parameters{5} = 0;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in radiobutton_no_inverted.
function radiobutton_no_inverted_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to radiobutton_no_inverted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_no_inverted
if get(hObject,'Value') == 1
    set(handles.radiobutton_yes_inverted,'Value',0)
    handles.parameters{5} = 0;
else
    set(handles.radiobutton_yes_inverted,'Value',1)
    handles.parameters{5} = 1;
end

% Update handles
guidata(hObject,handles)


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
        

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Delete figure
delete(handles.figure1)


% --- Executes on button press in pushbutton_default.
function pushbutton_default_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset the parameters
handles.parameters{1} = 300;
handles.parameters{2} = 100;
handles.parameters{3} = 1;
handles.parameters{4} = 0;
handles.parameters{5} = 0;

% Update the GUI
update_gui(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);
