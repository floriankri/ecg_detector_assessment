function varargout = rename_channel(varargin)
% RENAME_CHANNEL MATLAB code for rename_channel.fig
%      RENAME_CHANNEL by itself, creates a new RENAME_CHANNEL or raises the
%      existing singleton*.
%
%      H = RENAME_CHANNEL returns the handle to a new RENAME_CHANNEL or the handle to
%      the existing singleton*.
%
%      RENAME_CHANNEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RENAME_CHANNEL.M with the given input arguments.
%
%      RENAME_CHANNEL('Property','Value',...) creates a new RENAME_CHANNEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rename_channel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rename_channel_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help rename_channel

% Last Modified by GUIDE v2.5 06-May-2019 09:50:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rename_channel_OpeningFcn, ...
                   'gui_OutputFcn',  @rename_channel_OutputFcn, ...
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


% --- Executes just before rename_channel is made visible.
function rename_channel_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rename_channel (see VARARGIN)

% Adjust the current name
set(handles.text_current_name2,'string',varargin{2});

% Choose default command line output for rename_channel
handles.output = '';

% Get nice dialog locations
set(handles.figure_rename_channel,'Position',getnicedialoglocation(get(hObject,'Position'), get(hObject,'Units')));

% Set default button
handles.figure_rename_channel.setDefaultButton(handles.pushbutton_ok);

% Set the keypress function
set([handles.figure_rename_channel, handles.pushbutton_ok, handles.edit_new_name], 'KeyPressFcn', {@(hObject,eventdata)rename_channel('doKeypress',hObject,eventdata,guidata(hObject))});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rename_channel wait for user response (see UIRESUME)
uiwait(handles.figure_rename_channel);


% --- Outputs from this function are returned to the command line.
function varargout = rename_channel_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure_rename_channel);


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Define the output
handles.output = get(handles.edit_new_name,'string');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure_rename_channel);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure_rename_channel);


% --- Executes when user attempts to close figure_rename_channel.
function figure_rename_channel_CloseRequestFcn(hObject, eventdata, handles) %#ok
% hObject    handle to figure_rename_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press
function doKeypress(hObject, eventdata, handles) %#ok
% hObject    handle to figure_rename_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
switch eventdata.Key
    case 'escape'
        
        % Resume
        uiresume(handles.figure_rename_channel);
        
    case 'return'
        % Define the output
        handles.output = get(handles.edit_new_name,'string');

        % Update handles structure
        guidata(hObject, handles);

        % Resume
        uiresume(handles.figure_rename_channel);
        
    otherwise
        % Update handles structure
        guidata(hObject, handles);
end    


function edit_new_name_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to edit_new_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_new_name as text
%        str2double(get(hObject,'String')) returns contents of edit_new_name as a double

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_new_name_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to edit_new_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
