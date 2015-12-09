function varargout = overshooting_gui(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @overshooting_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @overshooting_gui_OutputFcn, ...
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


% --- Executes just before overshooting_gui is made visible.
function overshooting_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for overshooting_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

ybar_v = log(1800);
istar_v = log(0.025);
pstar_v = 0;
m_v = log(150);

theta_v = 0.06;
lambda_v = 0.80;
delta_v = 20;
sigma_v =  0.1;
alpha_v = 0.06;


%Set valores en GUI
ybar =findobj('Tag','ybar_i');
set(ybar,'String',exp(ybar_v));
istar =findobj('Tag','istar_i');
set(istar,'String',exp(istar_v));
pstar =findobj('Tag','pstar_i');
set(pstar,'String',exp(pstar_v));
m =findobj('Tag','m');
set(m,'String',exp(m_v));

theta =findobj('Tag','theta_i');
set(theta,'String',(theta_v));
lambda =findobj('Tag','lambda_i');
set(lambda,'String',(lambda_v));
delta =findobj('Tag','delta_i');
set(delta,'String',(delta_v));
alpha =findobj('Tag','alpha_i');
set(alpha,'String',(alpha_v));
sigma =findobj('Tag','sigma_i');
set(sigma,'String',(sigma_v));
% UIWAIT makes overshooting_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = overshooting_gui_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function alpha_i_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function value=getNumberValue(object)
    value =findobj('Tag',object);
    value = get(value,'String');
    value = eval(value);
    value;

%GRAFICA CURVAS DE FASE CON CONDICIONES BASE.
function pushbutton1_Callback(hObject, eventdata, handles)
ybar = log(getNumberValue('ybar_i'));
istar = log(getNumberValue('istar_i'));
pstar = log(getNumberValue('pstar_i'));
m = log(getNumberValue('m'));

theta = getNumberValue('theta_i');
lambda = getNumberValue('lambda_i');
delta = getNumberValue('delta_i');
sigma = getNumberValue('sigma_i');
alpha = getNumberValue('alpha_i');

limit = 6;
%Determinando la matriz A y b
A = [(-alpha)*(delta + sigma/lambda) alpha*delta; 1/lambda 0];
b = [alpha*((sigma*m)*(1/lambda) - (1 + (sigma*theta)*(1/lambda))*ybar) ; (1/lambda)*(theta*ybar - m) - istar];
Mss = -(A^-1)*b;

% Valores de Estado Estacionario
pss = Mss(1,1);
sss = Mss(2,1);

%Set valores en GUI
pss_o =findobj('Tag','pss');
set(pss_o,'String',pss);
sss_o =findobj('Tag','sss');
set(sss_o,'String',sss);

st = 0 : 1 : limit;

axes(handles.axes1);

% Ceroclina del sistema es pt = 
p_fase = (delta / ((sigma/lambda) + delta))*(st -sss) + pss;

plot([0 limit],[pss pss],'--b',st,p_fase,':r',sss,pss,'*m','LineWidth',2); hold on;

grid on;
title('Diagrama de fases para tipo de cambio y precios domesticos')
xlabel(handles.axes1,'Tipo de cambio')
ylabel(handles.axes1,'Precios domesticos')
legend('Curva de fase S','Curva de fase P','equilibrio de EE')


%GRAFICA EL RECORRIDO DADO S0 Y P0 COMO COND. INICIALES ELEGIDAS
function pushbutton2_Callback(hObject, eventdata, handles)
ybar = log(getNumberValue('ybar_i'));
istar = log(getNumberValue('istar_i'));
pstar = log(getNumberValue('pstar_i'));
m = log(getNumberValue('m'));

theta = getNumberValue('theta_i');
lambda = getNumberValue('lambda_i');
delta = getNumberValue('delta_i');
sigma = getNumberValue('sigma_i');
alpha = getNumberValue('alpha_i');

limit = 6;
%Determinando la matriz A y b
A = [(-alpha)*(delta + sigma/lambda) alpha*delta; 1/lambda 0];
b = [alpha*((sigma*m)*(1/lambda) - (1 + (sigma*theta)*(1/lambda))*ybar) ; (1/lambda)*(theta*ybar - m) - istar];

%Condiciones Iniciales
St = getNumberValue('s0');
Pt = getNumberValue('p0');

i=1;
while i<60;
    %agregamos el punto anterior
    if St < limit && Pt < limit && St >= 0 && Pt >= 0-2 
        
        %Plotting path
        axes(handles.axes1);
        plot(St,Pt,'*g');
        
        %Plotting cambios en S y P
        axes(handles.axes2);
        plot(i-1,St,'*r','LineWidth',1);
        grid on;
        xlabel(handles.axes2,'Time Step')
        ylabel(handles.axes2,'Tipo de cambio')
        
        axes(handles.axes3);
        plot(i-1,Pt,'*r','LineWidth',1);
        grid on;
        xlabel(handles.axes3,'Time Step')
        ylabel(handles.axes3,'Precios domesticos')

        smoothening = 0.05;
        %Calculamos pr?ximo punto partiendo de matrices
        Mt = (A*[Pt;St] + b)*smoothening;
        Pt = Mt(1,1) + Pt;
        St = Mt(2,1) + St;
        
        %pausamos el computo 0.15 segundos
        pause(0.15)
    else
    end
    i = i+1;
end;

%GRAFICA EL RECORRIDO POR SADDLE PATH HASTA EL ESTADO ESTACIONARIO DE LA RAIZ ESTABLE DADO S0.
function pushbutton4_Callback(hObject, eventdata, handles)
ybar = log(getNumberValue('ybar_i'));
istar = log(getNumberValue('istar_i'));
pstar = log(getNumberValue('pstar_i'));
m = log(getNumberValue('m'));

theta = getNumberValue('theta_i');
lambda = getNumberValue('lambda_i');
delta = getNumberValue('delta_i');
sigma = getNumberValue('sigma_i');
alpha = getNumberValue('alpha_i');

limit = 6;
%Determinando la matriz A y b
A = [(-alpha)*(delta + sigma/lambda) alpha*delta; 1/lambda 0];
b = [alpha*((sigma*m)*(1/lambda) - (1 + (sigma*theta)*(1/lambda))*ybar) ; (1/lambda)*(theta*ybar - m) - istar];
Mss = -(A^-1)*b;

% Valores de Estado Estacionario
Pss = Mss(1,1);
Sss = Mss(2,1);

%Calculando las ra?ces
trA = -alpha*(delta + (sigma/lambda));
DELTA = (alpha^2)*(delta + (sigma/lambda))^2 + 4*(alpha*delta/lambda);
tau2 = ((trA) - sqrt(DELTA))/2;

St = 0 : 1 : 4;
axes(handles.axes1);
% Brazo estable - saddle path
p_saddle = (St - Sss + (1/(lambda*tau2))*Pss)*(lambda*tau2);

plot(St,p_saddle,':y','LineWidth',3); hold on;
legend('Curva de fase S','Curva de fase P','equilibrio de EE','Brazo estable')

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
ybar = log(getNumberValue('ybar_i'));
istar = log(getNumberValue('istar_i'));
pstar = log(getNumberValue('pstar_i'));
m = log(getNumberValue('m'));

theta = getNumberValue('theta_i');
lambda = getNumberValue('lambda_i');
delta = getNumberValue('delta_i');
sigma = getNumberValue('sigma_i');
alpha = getNumberValue('alpha_i');

limit = 6;
%Condiciones Iniciales
St = getNumberValue('S0_saddle');

%Determinando la matriz A y b
A = [(-alpha)*(delta + sigma/lambda) alpha*delta; 1/lambda 0];
b = [alpha*((sigma*m)*(1/lambda) - (1 + (sigma*theta)*(1/lambda))*ybar) ; (1/lambda)*(theta*ybar - m) - istar];
Mss = -(A^-1)*b;

% Valores de Estado Estacionario
Pss = getNumberValue('pss')
Sss = getNumberValue('sss')

%Calculando las raices
trA = -alpha*(delta + (sigma/lambda));
DELTA = (alpha^2)*(delta + (sigma/lambda))^2 + 4*(alpha*delta/lambda);
tau2 = ((trA) - sqrt(DELTA))/2;

axes(handles.axes1);


% Condiciones iniciales Saddle Path
Pt =  (St - Sss)*(lambda*tau2) + Pss;

i=1;
while i<60;
    %agregamos el punto anterior
    if St < limit && Pt < limit && St >= 0 && Pt >= -2  
        
        %Plotting path
        axes(handles.axes1);
        plot(St,Pt,'*g');
        
        %Plotting cambios en S y P
        axes(handles.axes2);
        plot(i-1,St,'*r','LineWidth',1);
        grid on;
        xlabel(handles.axes2,'Time Step')
        ylabel(handles.axes2,'Tipo de cambio')
        
        axes(handles.axes3);
        plot(i-1,Pt,'*r','LineWidth',1);
        grid on;
        xlabel(handles.axes3,'Time Step')
        ylabel(handles.axes3,'Precios domesticos')

        smoothening = 0.1;
        %Calculamos proximo punto partiendo de matrices
        Mt = (A*[Pt;St] + b)*smoothening;
        Pt = Mt(1,1) + Pt;
        St = Mt(2,1) + St;
        
        %pausamos el computo 0.15 segundos
        pause(0.15)
    else
    end
    i = i+1;
end;

%SHOCKEA LA OFERTA MONETARIA Y GRAFICA EL RECORRIDO DE VUELTA AL ESTADO ESTACIONARIO.
function pushbutton5_Callback(hObject, eventdata, handles)
ybar = log(getNumberValue('ybar_i'));
istar = log(getNumberValue('istar_i'));
pstar = log(getNumberValue('pstar_i'));
m = log(getNumberValue('m_nuevo') + getNumberValue('m'));

%Shock a la oferta monetaria propagado a la GUI
m_o =findobj('Tag','m');
set(m_o,'String',exp(m));

theta = getNumberValue('theta_i');
lambda = getNumberValue('lambda_i');
delta = getNumberValue('delta_i');
sigma = getNumberValue('sigma_i');
alpha = getNumberValue('alpha_i');

limit = 6;
%Determinando la matriz A y b
A = [(-alpha)*(delta + sigma/lambda) alpha*delta; 1/lambda 0];
b = [alpha*((sigma*m)*(1/lambda) - (1 + (sigma*theta)*(1/lambda))*ybar) ; (1/lambda)*(theta*ybar - m) - istar];
Mss = -(A^-1)*b;

% Valores de Estado Estacionario
pss_old = getNumberValue('pss');
sss_old = getNumberValue('sss');
pss = Mss(1,1);
sss = Mss(2,1);

%Set valores en GUI
pss_o =findobj('Tag','pss');
set(pss_o,'String',pss);
sss_o =findobj('Tag','sss');
set(sss_o,'String',sss);

st = 0 : 1 : limit;

axes(handles.axes1);

% Ceroclina del sistema es pt = 
p_fase = (delta / ((sigma/lambda) + delta))*(st -sss) + pss;

plot([0 limit],[pss pss],'--b',st,p_fase,':r',sss,pss,'*m','LineWidth',2); hold on;

%Calculando las raices
trA = -alpha*(delta + (sigma/lambda));
DELTA = (alpha^2)*(delta + (sigma/lambda))^2 + 4*(alpha*delta/lambda);
tau2 = ((trA) - sqrt(DELTA))/2;

st = 0 : 1 : 5;
axes(handles.axes1);
% Brazo estable - saddle path
p_saddle = (st - sss + (1/(lambda*tau2))*pss)*(lambda*tau2);

plot(st,p_saddle,':y','LineWidth',3); hold on;

St = sss_old;
Pt = pss_old;


i=1;
while i<60;
    %Plotting cambios en S y P
    axes(handles.axes2);
    plot(i-1,St,'*r','LineWidth',1);
    grid on;
    xlabel(handles.axes2,'Time Step')
    ylabel(handles.axes2,'Tipo de cambio')

    axes(handles.axes3);
    plot(i-1,Pt,'*r','LineWidth',1);
    grid on;
    xlabel(handles.axes3,'Time Step')
    ylabel(handles.axes3,'Precios domesticos')
    
    %agregamos el punto anterior
    axes(handles.axes1);
    if i == 1
        s_first_shock = (sss - (1/(lambda*tau2))*pss) + (1/(lambda*tau2))*pss_old;
        %Plotting initial shock
        plot([sss_old s_first_shock],[pss_old pss_old],'g','LineWidth',4);
        St = s_first_shock;
        Pt = pss_old;
    else
        plot(St,Pt,'*g');
    end

    smoothening = 0.1;
    %Calculamos proximo punto partiendo de matrices
    Mt = (A*[Pt;St] + b)*smoothening;
    Pt = Mt(1,1) + Pt;
    St = Mt(2,1) + St;

    %pausamos el computo 0.15 segundos
    pause(0.15)

    i = i+1;
end;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla
axes(handles.axes2);
cla
axes(handles.axes3);
cla


function ybar_i_Callback(hObject, eventdata, handles)
% hObject    handle to ybar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ybar_i as text
%        str2double(get(hObject,'String')) returns contents of ybar_i as a double


% --- Executes during object creation, after setting all properties.
function ybar_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ybar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta_i_Callback(hObject, eventdata, handles)
% hObject    handle to theta_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta_i as text
%        str2double(get(hObject,'String')) returns contents of theta_i as a double


% --- Executes during object creation, after setting all properties.
function theta_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m_Callback(hObject, eventdata, handles)
% hObject    handle to m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m as text
%        str2double(get(hObject,'String')) returns contents of m as a double


% --- Executes during object creation, after setting all properties.
function m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function istar_i_Callback(hObject, eventdata, handles)
% hObject    handle to istar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of istar_i as text
%        str2double(get(hObject,'String')) returns contents of istar_i as a double


% --- Executes during object creation, after setting all properties.
function istar_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to istar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pstar_i_Callback(hObject, eventdata, handles)
% hObject    handle to pstar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pstar_i as text
%        str2double(get(hObject,'String')) returns contents of pstar_i as a double


% --- Executes during object creation, after setting all properties.
function pstar_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pstar_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_i_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_i as text
%        str2double(get(hObject,'String')) returns contents of sigma_i as a double


% --- Executes during object creation, after setting all properties.
function sigma_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delta_i_Callback(hObject, eventdata, handles)
% hObject    handle to delta_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_i as text
%        str2double(get(hObject,'String')) returns contents of delta_i as a double


% --- Executes during object creation, after setting all properties.
function delta_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_i_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_i as text
%        str2double(get(hObject,'String')) returns contents of lambda_i as a double


% --- Executes during object creation, after setting all properties.
function lambda_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function s0_Callback(hObject, eventdata, handles)
% hObject    handle to s0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s0 as text
%        str2double(get(hObject,'String')) returns contents of s0 as a double


% --- Executes during object creation, after setting all properties.
function s0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0 as text
%        str2double(get(hObject,'String')) returns contents of p0 as a double


% --- Executes during object creation, after setting all properties.
function p0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function S0_saddle_Callback(hObject, eventdata, handles)
% hObject    handle to S0_saddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of S0_saddle as text
%        str2double(get(hObject,'String')) returns contents of S0_saddle as a double


% --- Executes during object creation, after setting all properties.
function S0_saddle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S0_saddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m_nuevo_Callback(hObject, eventdata, handles)
% hObject    handle to m_nuevo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m_nuevo as text
%        str2double(get(hObject,'String')) returns contents of m_nuevo as a double


% --- Executes during object creation, after setting all properties.
function m_nuevo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m_nuevo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


