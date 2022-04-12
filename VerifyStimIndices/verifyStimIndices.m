function varargout = verifyStimIndices(varargin)
% VERIFYSTIMINDICES MATLAB code for verifyStimIndices.fig
%      VERIFYSTIMINDICES, by itself, creates a new VERIFYSTIMINDICES or raises the existing
%      singleton*.
%
%      H = VERIFYSTIMINDICES returns the handle to a new VERIFYSTIMINDICES or the handle to
%      the existing singleton*.
%
%      VERIFYSTIMINDICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VERIFYSTIMINDICES.M with the given input arguments.
%
%      VERIFYSTIMINDICES('Property','Value',...) creates a new VERIFYSTIMINDICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before verifyStimIndices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to verifyStimIndices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help verifyStimIndices

% Last Modified by GUIDE v2.5 25-Jun-2013 12:35:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @verifyStimIndices_OpeningFcn, ...
                   'gui_OutputFcn',  @verifyStimIndices_OutputFcn, ...
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


% --- Executes just before verifyStimIndices is made visible.
function verifyStimIndices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to verifyStimIndices (see VARARGIN)

% Choose default command line output for verifyStimIndices
handles.output = false;

handles.hB = varargin{1};
if ~(isa(handles.hB, 'BardFile') || isa(handles.hB, 'MacLabFile'))
    error('wrong input type')
end

set(handles.editIndices, 'String', num2str([handles.hB.IsStimCaptured handles.hB.StimIndices]))
set(handles.popupmenu1, 'String', {'none' handles.hB.ChName{:}})
update_axes(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes verifyStimIndices wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = verifyStimIndices_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1)


function editIndices_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to editIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIndices as text
%        str2double(get(hObject,'String')) returns contents of editIndices as a double


% --- Executes during object creation, after setting all properties.
function editIndices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushUpdate.
function pushUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_axes(handles)

% --- Executes on button press in pushDone.
function pushDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val  = str2num(get(handles.editIndices, 'String')); %#ok<ST2NM>
if isempty(val)
    iStim = [];
    isCaptured = [];
else
    iStim = val(:,2);
    isCaptured = val(:,1);
end
handles.hB.StimIndices = iStim;
handles.hB.IsStimCaptured = isCaptured;
handles.output = true;
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in pushCancel.
function pushCancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to pushCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = false;
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Updates the axes when called.
function update_axes(handles)
hB = handles.hB;

e = egm(hB, ':', hB.ChStim);
val  = str2num(get(handles.editIndices, 'String')); %#ok<ST2NM>
if isempty(val)
    iStim = [];
    isCaptured = [];
else
    iStim = val(:,2);
    isCaptured = val(:,1);
end

val = get(handles.popupmenu1, 'Value');
e2 = [];
if val ~= 1 %'non'
    val = val-1;
    ch2 = hB.ChName{val};
    e2 = egm(hB, ':', ch2);
end

xLim = get(handles.axes1, 'XLim');
yLim = get(handles.axes1, 'YLim');
u = get(handles.axes1, 'Userdata');

set(handles.axes1, 'NextPlot', 'replace')
plot(handles.axes1, e)
set(handles.axes1, 'NextPlot', 'add')

plot(handles.axes1, e2, 'g')

for i = 1:numel(iStim)
    if isCaptured(i)
        col = 'r';
    else
        col = 'k';
    end
    plot(handles.axes1, [iStim(i) iStim(i)], [-0.5 +0.5], 'Color',col , 'LineWidth' , 2)
end


if isfield(u,'hasRun')
    set(handles.axes1, 'XLim',xLim , 'YLim',yLim);
end
u.hasRun = true;
set(handles.axes1, 'Userdata',u )


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
update_axes(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
