function varargout = DVOR_Simulator(varargin)
%DVOR_SIMULATOR MATLAB code file for DVOR_Simulator.fig
%      DVOR_SIMULATOR, by itself, creates a new DVOR_SIMULATOR or raises the existing
%      singleton*.
%
%      H = DVOR_SIMULATOR returns the handle to a new DVOR_SIMULATOR or the handle to
%      the existing singleton*.
%
%      DVOR_SIMULATOR('Property','Value',...) creates a new DVOR_SIMULATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to DVOR_Simulator_OpineningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DVOR_SIMULATOR('CALLBACK') and DVOR_SIMULATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVOR_SIMULATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DVOR_Simulator

% Last Modified by GUIDE v2.5 24-Jan-2022 15:23:22

% Begin initialization code - DO NOT EDIT
%%


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DVOR_Simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @DVOR_Simulator_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before DVOR_Simulator is made visible.
function DVOR_Simulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for DVOR_Simulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DVOR_Simulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DVOR_Simulator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLmAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_Start.
function button_Start_Callback(hObject, eventdata, handles)
% hObject    handle to button_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%--------DVOR����������------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Ptable;

format long;
try
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
   ph_ref=get(handles.edit_AZ_Align,'String');  %30Hz���࣬���ڷ�λУ��
    
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
    az_align=str2double(ph_ref);
    
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
 
  fc=freq_value*1e6;
   fc1=fc;  %�ز�Ƶ�ʣ���λHz
  
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
 
   f_Lsb1=(fc1-9960);       %LSBƵ��
  f_Usb1=(fc1+9960); 
  
  wave_LSB=300/f_Lsb;
  wave_USB=300/f_Usb;
  k=2*pi/wave_L;   %��λ���������ڲ�����λ����
  k_LSB=2*pi/wave_LSB; 
  k_USB=2*pi/wave_USB;
  
  omega=30;            % �ز�����Ƶ��
  ph30=az_align;              % ��Ƶ��λ
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
%  fs=200e3;   %����Ƶ��,�������ԣ���4e6��6��8��10���Եõ��ԳƵ��ز�AM�ź� ����3��5e6��7��9���Գ�
 fs=4e6;
 
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%  noisy=randn(1,N+1);
 
  t=0:1/fs:T-1/fs;    %����2��ʱ��Ƭ���ڵ��ź�;
  csb_power=str2double(get(handles.edit_CSB_Power,'String'));     % CSB����,��λW
  Lsb_phase=str2double(get(handles.edit_LSB_Phase,'String'));     %LSB��λ����λ�ȡ�
  Lsb_Amp=str2double(get(handles.edit_LSB_AMP,'String'))/100;     %LSB����
  Usb_Amp=str2double(get(handles.edit_LSB_USB_Ratio,'String'))/100;   % USB����
  AM30=str2double(get(handles.edit_AM30,'String'))/100;      % 30Hz AM ���ƶ�
  CSB_A=sqrt(2*csb_power*R);      
  
  csb_sideband=AM30*sin(2*pi*omega*t+ph30*pi/180);
% clr_sideband_150=n*sin(2*pi*omega2*t);
   
%    FFT_30=sin(2*pi*omega*t);
   
csb_mod=csb_sideband;
  
  CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
  LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
  USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
  SB_signal=LSB+USB;
  sum_s=CSB+SB_signal;    %���ź�=�ز�+�ߴ�
  v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');   %��Ӹ�˹�������������źŵ�20dB
  
  else
   AM_signal=sum_s;
  end
  

  h=hilbert(AM_signal);   %����Hilbert�任������90��
%  yi=imag(h);       %Hilbert�任֮��õ��ĸ������鲿����Hilbert�任
%   xi=real(h);
 
  
  
  
 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 
%  am_env=AM_signal;
           %  reshape(z,1,10000);
%  am_env=zz;

 %%%%%%%%%%%%-----------------------����FFT���ܣ�--------------------------
%    amp1=mean(am_env);
%   dft1=am_env.*FFT_30;   %�ֶ����㣬30Hz��DFT�任��ԭ������30Hz�źŶԲ����źŽ��о��������Ӧ����������
%   amp2=sum(dft1);        %���ݹ�ʽ���㣬
%   amp3=amp2/(N/2);       % �õ�30HzƵ�׳ɷݵķ��ȣ���FFT��ȣ�һ�¡���֤��ϡ�
  
%  subplot(4,1,3);
%   plot(t,zzz);
%   ylabel('����');xlabel('ʱ��');title('����첨');
 

 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 %fs=1/(T/(N-1));  %����Ƶ��
%  ffss=fs*linspace(0,1,NFFT);
 
 ft=fft(am_env,NFFT);
 
 FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 %fh1=FH.*conj(FH);;%conj()������y�����Ĺ������ʵ���Ĺ������������
% ff=fs*(0:NFFT/2-1);% F F T �任���Ӧ��Ƶ�ʵ�����
 %b=FH(90);
 
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
  
  index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 
 
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);

  
  
  
  ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����


xxlable="ʱ��";        %get(handles.edit_XLable,'String');
yylable="����";         %get(handles.edit_YLable,'String');
  figure(3);
  
  subplot(3,2,1);
plot(t,CSB);
% xlabel(xxlable);2
ylabel(yylable);
title(ftitle);
  t_p=T;
 axis([0 t_p -CSB_A*(1+AM30) CSB_A*(1+AM30)]);
 
   subplot(3,2,2);
    ftitle="�ߴ�";
plot(t,SB_signal);
% xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -CSB_A CSB_A]);
   
   
   subplot(3,2,3);
    ftitle="�ռ�������ź�";
plot(t,AM_signal);
% xlabel(xxlable);
ylabel(yylable);
title(ftitle);
    t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);
   
    
   subplot(3,2,4);
    ftitle="AM�����  AM30="+num2str(AM30_MOD,'%1.4f')+"  AM9960="+num2str(AM9960_MOD,'%1.4f');
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);
  
  str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
 t_p=0.07;
 
 text(t_p/2,-30, str1,'Color','red','FontSize',8);
 str2={num2str(FH(1)),  num2str(FH(index_30),'%1.4f'),  num2str(FH(index_9960),'%1.4f')};
  text(t_p/2+40*t_p/80,-30,str2,'Color','red','FontSize',8);
 
   subplot(3,2,5);
    ftitle="Ƶ��ͼ";
  freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
  YY=fftshift(FH);
  plot(freqaxis,YY);

xlabel("Ƶ��");
ylabel("����");
title(ftitle);
 
%  axis([-N N -0.2*CSB_A 2*CSB_A]);
grid on

%%%%%%%%%%------------������Ϻ��� ����-------------------
blending_f=get(handles.popupmenu_BlendingFunction,'Value');
w = 2*T_on;     %1/1440*2=1/720  �ߴ���ʱ�䣬Ϊ��������
t_sb=fix(w/T*length(t)); %�ߴ��򿪵�ʱ��Ƭ�����������ݳ��ȣ�
if (mod(t_sb,2))~=0     %ȡż��
    t_sb=t_sb+1;
end

shift_t=t_sb/2;    %ƽ��50%�����
shift_sb=zeros(1,shift_t); %����һ��ƽ�������,�����桢ż��

t_csb=t;

% t=t(1:t_sb);





switch blending_f
    case 1   % 'COS^2/SIN^2'
      b_cos=cos(2*pi*(1/(4*T_on))*t).*cos(2*pi*(1/(4*T_on))*t);  %��Ϻ�����������T_on�������� T_on�Ǳߴ���ʱ��
      b_sin=sin(2*pi*(1/(4*T_on))*t).*sin(2*pi*(1/(4*T_on))*t); 
    case 2    %'COS/SIN'
        b_cos=cos(2*pi*(1/(4*T_on))*t);
         b_sin=sin(2*pi*(1/(4*T_on))*t);
%     case 'CSC'
        
    case 3   % 'SQUARE'
        b_cos=1;
        b_sin=1;
end




% b=cos(2*pi*720*t).*cos(2*pi*720*t);
% b=1;
% % % % % x_odd = rectpuls(t,w);   %�����������壬���Ϊw
% % % % % x_even=rectpuls(t+T_on/2,w);
% % % % % 
% % % % % x2_odd=x_odd.*b_cos;              % ������ϵ��ư���
% % % % % x2_even=x_even.*b_sin;


%%
%
% x3=csb.*x2;
% 
% tpast = -T_on; %1/1440; ��ǰһ��T_on
% xpast = rectpuls(t-tpast,w).*b_cos;
% 
% tfutr=T_on;  %1/1440;   �ͺ�һ��T_on
% xfutr = rectpuls(t-3*tfutr,w).*b_cos;
% plot(t,x,t,xpast,t,xfutr,t,x2,t,x3);

% y=square(2*pi*30*t,T_on/T*100);
% [a,b]=size(y);
% for i=1:b
%     if y(i)<0
%         y(i)=0;
%     end
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
%   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
% CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);

CSB_func=2*pi*fc1*t_csb;

LSB_func=2*pi*f_Lsb1*t;
USB_func=2*pi*f_Usb1*t;


 % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




% % % % % % % % % % % % % %   antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
% % % % % % % % % % % % % %     couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
% % % % % % % % % % % % % %     counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
% % % % % % % % % % % % % %   reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
% % % % % % % % % % % % % %   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
% % % % % % % % % % % % % %    CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
% % % % % % % % % % % % % %        SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
% % % % % % % % % % % % % %      FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
% % % % % % % % % % % % % %      start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %      end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %        start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ����
% % % % % % % % % % % % % %        end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %��ʼ����
% % % % % % % % % % % % % %    FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
% % % % % % % % % % % % % %     simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�

%%%%%%%%%%%%%%%%%%%%%%Բ�ܷ���%%%%%%%%%%%%%
start_angle=start_range;
stop_angle=end_range;
hwait=waitbar(0,'processing...0%','name','�ߴ��ź������У�please wait>>>');


sbs_total=T/T_on;
steps=sbs_total/100;
CSB_X=0;
CSB_Y=0;
CSB_Z=CSB_H;
CSB_elev=atan(CSB_H/counterpoint_R)*180/pi;   %CSB������Ե��������ǣ���Ϊ�ж����Է��������Ǵ����Ϊ����ο�

% LSB_ANT=zeros(sbants,length(t));
% LSB_ANT_IMG=zeros(sbants,length(t));
% USB_ANT=zeros(sbants,length(t));
% USB_ANT_IMG=zeros(sbants,length(t));

% LSB_ANT=zeros(1,t_sb);
% LSB_ANT_IMG=zeros(1,t_sb);
% USB_ANT=zeros(1,t_sb);
% USB_ANT_IMG=zeros(1,t_sb);

EVEN_ANT=zeros(1,t_sb);
EVEN_LSB=zeros(1,t_sb);
EVEN_USB=zeros(1,t_sb);

ODD_ANT=zeros(1,t_sb);
ODD_LSB=zeros(1,t_sb);
ODD_USB=zeros(1,t_sb);




sb_low=get(handles.uitable1,'Data');
sb_high=get(handles.uitable2,'Data');
sb_all=cat(1,sb_low,sb_high);
[a,b]=size(sb_all);
sb_basic_data=zeros(sbants,6);   %%1:X,2:Y,3:Z,4:A,5:P,6.ON/OFF  �ߴ����߻����������ݾ���

sb_basic_data_image_1_counterpoint=zeros(sbants,6);    %�Է������Ϊ��׼�ľ���
sb_basic_data_image_2_ground=zeros(sbants,6);          %�Դ��Ϊ��׼�ľ���

if ~isempty(sb_all{1,1})
for b=1:sbants     %���ɱߴ���������  1�����ߺţ�2���Ƕȣ�3�����룬4���߶ȣ�5�����ȣ�6����λ��7������
    sb_z=sb_all{b,4};
    sb_r=sb_all{b,3};
    sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}/100;
    sb_p=sb_all{b,6};
    sb_onoff=sb_all{b,7};
    
    sb_x=sb_r*sin(-sb_ang*pi/180);
    sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4};
%     sb_r=sb_all{b,3};     %ֻ��Z��ı䣬����X��Y������һ���ġ�
%     sb_ang=sb_all{b,2};
     sb_a=sb_all{b,5}*reflect_rate_couterpoint/100;
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data_image_1_counterpoint(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4}-2*couterpoint_H;
%     sb_r=sb_all{b,3};
%     sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}*reflect_rate_ground/100;   %���Ƿ���ϵ��
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
  sb_basic_data_image_2_ground(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
end
else
    set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',"���ȳ�ʼ���ߴ����������");
    return;
    
end

%%

%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Բ�ܷ��з���  ������������������������������
%%
mode_sel=get(handles.radiobutton_Circle,'Value');  %Բ�ܷ���ѡ��ť
if mode_sel==1       %ģʽѡ��=1����ʾԲ�ܷ���
simulation_range=start_angle:simulate_step:stop_angle;
step_c=length(simulation_range);

results=zeros(step_c,6);


for sim_step=1:step_c     %��Բ�ܵĲ�����ʼ���棬ÿ�����������0.1��ʱ���ڽ��յ����ź� 
    tic;
    ang=start_angle+(sim_step-1)*simulate_step;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=start_H-couterpoint_H;
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
% %     fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
    fly_d=FlySimulate_Circle;      %%  ����ϵԭ���ڵ������ġ�Y��ָ��1������������X��ָ�򶫷�
    fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
    
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
%     tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));   %%���㣨fly_x,fly_y,fly_z)�뾵������(0,0,-CSB_Z)���ߵĽ����Ƿ��ڵ��������ж�
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end

%                reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
%   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                  reflect_rate=reflect_rate_couterpoint;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                  
                   ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
                  
              end
   
     
    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+delta_csb_p;

  loss_csb=32.45+20*log10(fc)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
  CSB= csb_Loss*CSB_A*(1+csb_mod).*cos(CSB_func+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
    CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*cos(CSB_func+csb_img_p); %*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
   end
  
    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�д�����
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
     %receiver_SB�ӳ�����8696��
     %Բ�ܷ��С�
    [ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb,f_Usb,Lsb_phase,handles,counterpoint_R,hwait);
  
  
%%  %%  


%  for i=1:sbs_total   %    sbs_total=T/T_on  �ܵķ���ߴ��źţ��ڹ涨��ʱ��T�ڣ����� �ıߴ�������ÿ48��50һ��ѭ��
%      
% %     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
% %     y=SB_ANT(i,:);
% %     [a,b_cos]=size(y);
% % for k=1:b_cos
% %     if y(k)<0
% %         y(k)=0;
% %     end
% % end
% %   


%      % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
% %   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
%     
%     if mod(i,2)==0      %ż������
%         
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w); %���߿�ͨʱ�䣺һ����������
%          for pp=1:length(b_even)-1
%              if pp==1 && b_even(pp)==1
%                  even_start=1;
%              else
%              if b_even(pp)==0 && b_even(pp+1)==1    %�ҳ�ÿ�����߷������ֹʱ��
%                  even_start=pp+1;
%              end
%              end
%              
%              if pp==length(b_even)-1 && b_even(pp+1)==1
%                  even_stop=pp+1;
%              else
%                 if b_even(pp)==1 && b_even(pp+1)==0
%                  even_stop=pp;
%                 end
%              end
%              
%          end
%            
%              
% 
%          b_even_sb=b_even(even_start:even_stop).*b_sin(even_start:even_stop);  %��Ϻ���
%           
%          
% % b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w).*b_sin*Lsb_Amp*CSB_A; %���߿�ͨʱ�䣺һ����������  
% 
%          
%       sbnum_LSB=mod(i,sbants);
%       if sbnum_LSB==0
%           sbnum_LSB=sbants;    %�ߴ������1��ʼ������48��MOD����Ϊ0����ȡ48
%       end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
%                     
% %                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*cos(LSB_func(even_start:even_stop)+phase_error);
% 
%                     
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                         EVEN_LSB=zeros(1,even_stop-even_start);
%                            
% %                         
%                           
%                        else
%                            EVEN_LSB=LSB_f.*b_even_sb;
%                            
%                              
% %                          
%                        end
%                     %%%%%%%%      �ϰ��ﷴ����棬         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R   %sb_elev_max_lsb
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                   
%                   
%                   reflect_rate=reflect_rate_couterpoint;
%              
%                   
%               
%                   
%                else
% %                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                        sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                       sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                       sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                       sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                       
%                       reflect_rate=reflect_rate_ground;
%                       
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
%                 end
%                
%                
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*reflect_rate*sb_a*cos(LSB_func(even_start:even_stop)+phase_error);
% 
%                       %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                          if checkSB_OBS   %���ڵ�
%                              
%                              EVEN_LSB=EVEN_LSB+zeros(1,even_stop-even_start);
% %                         
%                           
%                          else
%                               EVEN_LSB=EVEN_LSB+LSB_f.*b_even_sb;
%                              
% %                          
%                         end
%                         
%                        
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                     
%                 
%                     
%       else  %��� ���߲����ã���������档
%                             
%           EVEN_LSB=zeros(1,even_stop-even_start);
%         
%           
%       end     %����������ã�LSB��IF���
%       
%     
%                      
%        %%%%%%%%%  ������USB    ************************************************************************      
%       
%       if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%        end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,������USB����������      
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*cos(USB_func(even_start:even_stop)+phase_error);                  
%                             
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                       if checkSB_OBS   %���ڵ�
%                                
%                              EVEN_USB=zeros(1,even_stop-even_start);
% %                         
%                           
%                          else
%                               EVEN_USB=USB_f.*b_even_sb;
%                              
% %                          
%                        end
%                         
%                         
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                    
%               
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%             
%                     
%                 if rr<= counterpoint_R  %sb_elev_max_usb
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                           reflect_rate=reflect_rate_couterpoint;
%                           
%                else
% %                    if  fly_angle<sb_elev_min_usb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                             reflect_rate=reflect_rate_ground;
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
% %                        
% %                    end
%                 end     %����ѡ���IF���
%                 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A* reflect_rate*sb_a*cos(USB_func(even_start:even_stop)+phase_error); 
%                    
%                              
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        
%                        if checkSB_OBS   %���ڵ�
%                                 
%                              EVEN_USB=EVEN_USB+zeros(1,even_stop-even_start);
% %                         
%                           
%                          else
%                               EVEN_USB=EVEN_USB+USB_f.*b_even_sb;
%                              
% %                          
%                         end 
%                         
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      
%             
%       else   %��� ���߲����ã�������
%              EVEN_USB=zeros(1,even_stop-even_start);
%           
%           
%            
%              
%       end    %�����������
%       
%            if i==2
%                EVEN_ANT=EVEN_LSB+EVEN_USB;
%            else
%                EVEN_ANT=[EVEN_ANT EVEN_LSB+EVEN_USB];
%            end    
%           
% 
%       
%     
%     else                    %%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         
%         
% %        b_odd= rectpuls(t,w);
%          b_odd= rectpuls(t-fix(i/2)*w,w);
%          
%          
%           for pp=1:length(b_odd)-1
%               if pp==1 && b_odd(pp)==1
%                   odd_start=1;
%               else
%                   
%              if b_odd(pp)==0 && b_odd(pp+1)==1    %�ҳ�ÿ�����߷������ֹʱ��
%                  odd_start=pp+1;
%              end
%               end
%               
%               if pp==length(b_odd)-1 && b_odd(pp+1)==1
%                   odd_stop=pp+1;
%               else
%                if b_odd(pp)==1 && b_odd(pp+1)==0  
%                  odd_stop=pp;
%                end
%               end
%              
%          end
%          b_odd_sb=b_odd(odd_start:odd_stop).*b_cos(odd_start:odd_stop);  %��Ϻ���
%          
%          
%          
%       sbnum_LSB=mod(i,sbants);
%          if sbnum_LSB==0
%           sbnum_LSB=sbants;
%          end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*cos(LSB_func(odd_start:odd_stop)+phase_error);
% 
%                               
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                         
%                         if checkSB_OBS   %���ڵ�
%                              
%                          ODD_LSB=zeros(1,odd_stop-odd_start);
%                            
% %                         
%                           
%                        else
%                            ODD_LSB=LSB_f.*b_odd_sb;
%                              
% %                          
%                         end 
%                         
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                          tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R   %sb_elev_max_lsb
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                    reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                    reflect_rate=reflect_rate_ground	;
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
%                 end
%                       
% %                if fly_angle>sb_elev_max_lsb
% %                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                else
% %                    if  fly_angle<sb_elev_min_lsb
% %                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
% %                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
% %                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
% %                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
% %                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
% %                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*reflect_rate*CSB_A*sb_a*cos(LSB_func(odd_start:odd_stop)+phase_error);
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                         
%                      if checkSB_OBS   %���ڵ�
%                          
%                                ODD_LSB=ODD_LSB+zeros(1,odd_stop-odd_start);
% %                         
%                           
%                          else
%                               ODD_LSB=ODD_LSB+LSB_f.*b_odd_sb;
%                              
%                              
%                              
% %                          
%                       end 
%                         
%                         
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                     
%                    
%                     
%       else
%            ODD_LSB=zeros(1,odd_stop-odd_start);
%           
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ���� ��USB���侵��**********************************************     
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*cos(USB_func(odd_start:odd_stop)+phase_error); %*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                    
%                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                         
%                         
%                      if checkSB_OBS   %���ڵ�
%                              
%                                 ODD_USB=zeros(1,odd_stop-odd_start);
%                            
% %                         
%                           
%                        else
%                            ODD_USB=USB_f.*b_odd_sb;
%                              
% %                          
%                       end 
%                        
%                        
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%              
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                      
% %                     if sbnum_LSB<=sbants/2
% %                     sbnum_USB=sbnum_LSB+sbants/2;
% %                     else
% %                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
% %                     end
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R   %sb_elev_max_lsb
%                 
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                           reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_usb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                       reflect_rate=reflect_rate_ground;
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
% %                        
% %                    end
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*reflect_rate*CSB_A*sb_a*cos(USB_func(odd_start:odd_stop)+phase_error); %*cos(phase_error)-sin(USB_func)*sin(phase_error));
%               %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                      
%                         if checkSB_OBS   %���ڵ�
%                             
%                             ODD_USB=ODD_USB+zeros(1,odd_stop-odd_start);
% %                         
%                           
%                          else
%                               ODD_USB=ODD_USB+USB_f.*b_odd_sb;
% %                          
%                        end 
%                         
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                    
%                   
%       else
%              ODD_USB=zeros(1,odd_stop-odd_start);
%          
%       end    %�����������
%        
% %          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
% %          USB_ANT(i,:)=USB.*b_odd.*b_cos;
%     if i==1
%         ODD_ANT=ODD_LSB+ODD_USB;
%     else
%         ODD_ANT=[ODD_ANT ODD_LSB+ODD_USB];
%     end
% 
%  
%     
%     
%     end
%     
%     
%     
%     PerStr=fix(i/steps);
%     waitstr=['processing.....',num2str(PerStr),'%','�� ',num2str(sbnum_LSB),' �űߴ�����'];
%     waitbar(i/sbs_total,hwait,waitstr);
% %     pause(0.0005);
%  end   
%%
 
 

 
%% %%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
%%   
  %%   %%%%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
  a=length(EVEN_ANT);
  b=length(ODD_ANT);
  if a>b
      ODD_ANT=[ODD_ANT zeros(1,a-b)];
  else
      EVEN_ANT=[EVEN_ANT zeros(1,b-a)];
  end
   sum_s=CSB+CSB_IMG+EVEN_ANT+ODD_ANT; 
%    sum_s=CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG); %um(USB_ANT)+sum(USB_ANT_IMG);


   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
%  am_env=AM_signal;
 


 
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
 index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);
 ang30=angle(ft(index_30))*180/pi;
 
 
 
   y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,���ݷ��ȼ��㹦�ʣ���1/2 A^2/R,��λ��dBW,����30������dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
   format long;
   %%
%                            N=new_nfft;
% %����Ƶ��
% fs=rtlsdr_fs;
% 
% %�����˲���������Ƶ��
% fc_lpf=30;
% fc_hpf=8000;
% fp_bandpass=[8000 10000];
% fc_stop=[200 400];
% 
% %�Բ���Ƶ�ʵ�һ�㣬��Ƶ�ʹ�һ��
% wn_lpf=fc_lpf*2/fs;
% wn_hpf=fc_hpf*2/fs;
% wn_bandpass=fp_bandpass*2/fs;
% wn_stop=fc_stop*2/fs;
% 
% %����fir1�������FIR�˲���
% b_lpf=fir1(N-1,wn_lpf);
% b_hpf=fir1(N-1,wn_hpf,'high');
% b_bandpass=fir1(N-1,wn_bandpass,'bandpass');
% b_stop=fir1(N-1,wn_stop,'stop');
% 
% yy=filter(b_lpf,1,yy_data);
%                                 
%  DC=mean(yy);
%  a30=yy-DC;
%  
%%                                
                    rtl_fft  =ft;
                    rtlsdr_data_fft=FH;
                      new_nfft=NFFT  ;        
                     rtlsdr_fs=fs;           
                      LOoffset=0;             
                                     span=30;
                                    [A30,index30]=max(rtlsdr_data_fft(maxId+(30-span)/freq_rev+1:maxId+(30+span)/freq_rev));

                                          
                                      ph30AM= angle(rtl_fft(maxId+(30-span)/freq_rev+index30))*180/pi; 
 
                                    fftDC=rtl_fft;
                                    fftDC(maxId+1:length(fftDC))=0;
                                    ss=ifft(fftDC,new_nfft);
                                   sss=real(ss);
%                                    [pp,locs]=findpeaks(sss);
                                   peaks=mean(sss);
%                                    peaks=median(pp);
%                                    peaks=max(pp);
                                   DC1=peaks;
                                                                     
                                 
                                   fft30= rtl_fft;
                                    span=30;
                                      fft30(maxId)=0;
                                     fft30(maxId:(30-span)/freq_rev)=0;
                                     fft30((30+span)/freq_rev+1:length(fft30)-(30+span)/freq_rev)=0;

                                     ss=ifft(fft30,new_nfft);

                                     am30_env=real(ss);
                                     
%                                      [pks,locs] = findpeaks(avSpots,'MinPeakDistance',6);
%                                        findpeaks(select,Fs,'MinPeakHeight',1)
                                     [pp,locs]=findpeaks(am30_env);
%                                    peaks=mean(pp);
                                   peaks=median(pp);
%                                    peaks=max(pp);
%                                    peak_int=trapz(am30_env(locs),pp);
                                   
                                   AM30=peaks/DC1;%FindAmp(9960,nfft,sss);
                                     

                                   
                                           
                                    fft1020= rtl_fft;
                                    span=300;
                                     fft1020(maxId:(1020-span)/freq_rev)=0;
                                     fft1020((1020+span)/freq_rev+1:length(fft1020)-(1020+span)/freq_rev)=0;
                                     fft1020(length(fft1020)-(1020-span)/freq_rev:end)=0;
                                   
                                     ss=ifft(fft1020,new_nfft);

                                     am1020_env=real(ss);
                                     [pp,locs]=findpeaks(am1020_env);
%                                    peaks=mean(pp);
                                   peaks=median(pp);
%                                    peaks=max(pp);
                                   
                                   AM1020=peaks/DC1;%FindAmp(9960,nfft,sss);
                                   
                                   
                                   

                                     span=660;
                                     fft9960=rtl_fft;
                                     

                                     fft9960(maxId:(9960-span)/freq_rev)=0;
                                     fft9960((9960+span)/freq_rev+1:length(fft9960)-(9960+span)/freq_rev)=0;
                                     fft9960(length(fft9960)-(9960-span)/freq_rev:end)=0;
                                     
                                     

                                      ss=ifft(fft9960,new_nfft);
                                   am9960_env=real(ss);
                             

                                   span=2000;
                                   
                            [max_f,min_f,mean_f,max_indix,min_indix]=frequencycounter(am9960_env(span:length(am9960_env)-span),rtlsdr_fs );
                             max_f=max_f*(1-LOoffset/rtlsdr_fs);
                             min_f=min_f*(1-LOoffset/rtlsdr_fs);
                             
                            
                             [pp,locs]=findpeaks(am9960_env(span:length(am9960_env)-span));
%                                    peaks=mean(pp);
                                   peaks=median(pp);
%                                    peaks=max(pp);
                                   
                                   AM9960=peaks/DC1;%FindAmp(9960,nfft,sss);
                            
                                   

                            
%%  FM Demodulation
fmdemod=am9960_env;
mmm=length(fmdemod);
% fmdemod(fmdemod>peaks/2)=peaks/2;
% fmdemod(fmdemod<-peaks/2)=-peaks/2;
t=0:1/rtlsdr_fs:mmm/rtlsdr_fs-1/rtlsdr_fs;
phase=angle(hilbert(fmdemod).*exp(-1i*2*pi*9960*t));
phi=unwrap(phase);
dem=(1/(2*pi*16))*diff(phi)/(1/rtlsdr_fs);
% dem(length(t))=dem(length(t)-1);
% figure(33);
% plot(t,dem);
%% --------------------------------
 
 
                                    
 
                                     fm30_env=real(dem);  %-mean(real(sss));
                                     
                                     
               
                                     
                            


                                     fm30FFT=fft(fm30_env);

                                   fm30FFTAMP=abs(fm30FFT);
                                   span=30;
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-span)/freq_rev+1:maxId+(30+span)/freq_rev));
                                    

                                   id30=id30+maxId+(30-span)/freq_rev;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;

                                      deg=ph30FM-ph30AM;
                                    
    

                                   if deg<0
                                        deg=deg+360;
                                    end
                                if deg>360
                                        deg=deg-360;
                                end
                                    
 
                                      fmi=(max_f- min_f)/2/30;
                                     
                                      

                                    
                                    az_error=deg-ang;  %���㷽λ���
                                    
                                    az_error=az_error-180; 
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
%                                        
                                      

                                          
                                          VOR_AZ=num2str(az_error);
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     results(sim_step,:)=[ang,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����
% xxlable="ʱ��";        %get(handles.edit_XLable,'String');
% yylable="����";         %get(handles.edit_YLable,'String');
%   figure(3);
  
%  
%    
%    subplot(3,2,3);
%     ftitle="�ռ�������ź�";
% plot(t,AM_signal);
% % xlabel(xxlable);
% ylabel(yylable);
% title(ftitle);
%     t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
   
    figure(3);
   subplot(3,2,6);
%     ftitle="AM������-AM30MOD="+VOR_30HzAM+"  AM9960MOD="+ VOR_9960HzAM;
      ftitle="AZ: "+num2str(ang)+" | AM30="+VOR_30HzAM+" | AM9960="+ VOR_9960HzAM+" | AZ_E_r_r="+VOR_AZ+"|FMI="+VOR_FMI;
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle,'FontSize' ,9);
%    t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
  
%   str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
%  t_p=0.07;
%  
%  text(t_p/2,1, str1,'Color','red','FontSize',8);
%  str2={num2str(FH(1)),  num2str(AM30_MOD,'%1.4f'),  num2str(AM9960_MOD,'%1.4f')};
%   text(t_p/2+33*t_p/200,1,str2,'Color','red','FontSize',8);
% 
% 
%    figure(10);
%     ftitle="Ƶ��ͼ";
%   freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
%   YY=fftshift(FH);
%   plot(freqaxis,YY);
% 
% xlabel("Ƶ��");
% ylabel("����");
% title(ftitle);
% grid on
%  

 
 
 
 
 %%
 
 
 
%  waitbar(0,hwait,'0%');
time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
    hwait.Name=['Left: ',num2str(fix(time_left)),'s ', num2str(time_over),'s/��',num2str(PerStr),'%'];
% pause(0.00005);
 
end    %����ѭ��������%
 close(hwait);
 
%  toc
 
% LSB_ANT1=y.*LSB;
% figure(4);
% for k=1:sbs_total
%     plot(t,USB_ANT(k,:));
% hold on;
% end
% plot(t,b_sin);
% plot(t,b_cos);
% 
% 
%  t_p=T;
%  axis([-t_p t_p -2*CSB_A 2*CSB_A]);
% grid on;

 %%%%% results(sim_step,:)=[ang,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%����������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
     
 figure(11);
 cc=results(:,1);
 subplot(5,1,1);
  plot(cc,results(:,3));
     ftitle="VOR��λ��� Vs. �Ƕ�";
xlabel("Circle");
ylabel("AZ error");
title(ftitle);

 subplot(5,1,2);
  plot(cc,results(:,4));
     ftitle="30HzAM Vs. �Ƕ�";
xlabel("Circle");
ylabel("30HzAM");
title(ftitle);

 subplot(5,1,3);
  plot(cc,results(:,5));
     ftitle="9960HzAM Vs. �Ƕ�";
xlabel("Circle");
ylabel("9960HzAM");
title(ftitle);

 subplot(5,1,4);
  plot(cc,results(:,2));
     ftitle="RF LEVEL Vs. �Ƕ�";
xlabel("Բ��");
ylabel("RF LEVEL��dBm)");
title(ftitle);


 subplot(5,1,5);
  plot(cc,results(:,6));
     ftitle="��Ƶָ��FMI Vs. �Ƕ�";
xlabel("Բ��");
ylabel("FMI");
title(ftitle);

end



 % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




% % % % % % % % % % % % % %   antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
% % % % % % % % % % % % % %     couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
% % % % % % % % % % % % % %     counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
% % % % % % % % % % % % % %   reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
% % % % % % % % % % % % % %   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
% % % % % % % % % % % % % %    CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
% % % % % % % % % % % % % %        SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
% % % % % % % % % % % % % %      FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
% % % % % % % % % % % % % %      start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %      end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %        start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ����
% % % % % % % % % % % % % %        end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %��ʼ����
% % % % % % % % % % % % % %    FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
% % % % % % % % % % % % % %     simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
%������ʼ�����ά����
start_z=start_H-couterpoint_H;
start_d=sqrt(start_range^2+start_z^2);
start_x=start_d*sin(FlySimulate_Radial*pi/180);
start_y=start_d*cos(FlySimulate_Radial*pi/180);

stop_z=end_H-couterpoint_H;
stop_d=sqrt(end_range^2+stop_z^2);
stop_x=stop_d*sin(FlySimulate_Radial*pi/180);
stop_y=stop_d*cos(FlySimulate_Radial*pi/180);

distan=sqrt((stop_z-start_z)^2+(stop_x-start_x)^2+(stop_y-start_y)^2);
dirVector=[stop_x-start_x,stop_y-start_y,stop_z-start_z]/distan;

step_select=get(handles.checkbox_D_H,'Value');  %�����Ļ�׼���Ǹ߶Ȼ��Ǿ���


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    �������   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_sel=get(handles.radiobutton_Radial,'Value');  %����ѡ��ť�����ѡ�У������1

if mode_sel==1   %ѡ�������ģʽ��
            simulation_range=0:simulate_step:distan;
step_c=length(simulation_range);

results=zeros(step_c,6);


for sim_step=1:step_c
    tic;
    new_d=(sim_step-1)*simulate_step; %�������� 
    
    newpoint=[start_x,start_y,start_z]+new_d.*dirVector; %�µ�����
    ang=FlySimulate_Radial;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=newpoint(3);
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(newpoint(1)^2+newpoint(2)^2+newpoint(3)^2);   %  FlySimulate_Circle^2-fly_z^2);
    fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%               tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    
      rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end
               if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                  reflect_rate=reflect_rate_couterpoint;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                 
                   ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
                  
              end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ ��2210�С�

    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+ delta_csb_p;

  loss_csb=32.45+20*log10(fc1)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc1)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
              %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
  CSB= csb_Loss*CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
     CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*(cos(CSB_func)*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
   end
  
    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�д�����
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
                 
 
   %  receiver_SB  function ��8696��
   %   ������з��� 
   [ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb,f_Usb,Lsb_phase,handles,counterpoint_R,hwait);
   
    
  %%  


%  for i=1:sbs_total   %sbants*30  �ܵķ���ߴ��ź�����
%      
% %     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
% %     y=SB_ANT(i,:);
% %     [a,b_cos]=size(y);
% % for k=1:b_cos
% %     if y(k)<0
% %         y(k)=0;
% %     end
% % end
% %   
%      % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
% %   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
%     
%     if mod(i,2)==0      %ż������
%         
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%          
%       sbnum_LSB=mod(i,sbants);
%       if sbnum_LSB==0
%           sbnum_LSB=sbants;
%       end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                       %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                             LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%           
%               if rr<=counterpoint_R
%                 
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                     reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                     reflect_rate=reflect_rate_ground ;
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*reflect_rate*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                            LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT_IMG(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                    
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,������USB����������      
%                     
%                      sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%              
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                                 
%                           sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                             reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_usb
%                           sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                             reflect_rate=reflect_rate_ground;
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
% %                        
% %                    end
%                end     %����ѡ���IF���
%                 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*reflect_rate*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT_IMG(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%               
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%       
%               
%                
%           
% 
%     
%     else                    %%%%%%%%      ------��������
%           
%        b_odd= rectpuls(t-fix(i/2)*w,w);
%        
%         
%       sbnum_LSB=mod(i,sbants);
%          if sbnum_LSB==0
%           sbnum_LSB=sbants;
%          end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                 LSB_ANT(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                     reflect_rate=reflect_rate_couterpoint;
%               else
%                  
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                     reflect_rate=reflect_rate_ground;
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% % %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*reflect_rate*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT_IMG(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ���� ��USB���侵��     
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%              
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                      
% %                     if sbnum_LSB<=sbants/2
% %                     sbnum_USB=sbnum_LSB+sbants/2;
% %                     else
% %                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
% %                     end
%                      tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                 
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                             reflect_rate=reflect_rate_couterpoint;
%                           
%               else
%                   
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                             reflect_rate=reflect_rate_ground;
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
% %                        
% %                    end
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                    
%                        phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                     
%                     
%                     
%                     
%                     
%                     
%                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                                USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                                  USB_ANT_IMG(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                 
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%        
% %          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
% %          USB_ANT(i,:)=USB.*b_odd.*b_cos;
%     
%       
% 
%     
%     
%     
%     end
%     
%     PerStr=fix(i/steps);
%     waitstr=['processing.....',num2str(PerStr),'%'];
%     waitbar(i/sbs_total,hwait,waitstr);
% %     pause(0.0005);
%  end
 %%
 
 
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
%%   
  %%   %%%%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
   a=length(EVEN_ANT);
  b=length(ODD_ANT);
  if a>b
      ODD_ANT=[ODD_ANT zeros(1,a-b)];
  else
      EVEN_ANT=[EVEN_ANT zeros(1,b-a)];
  end
   sum_s=CSB+CSB_IMG+EVEN_ANT+ODD_ANT; 
%    sum_s=CSB+CSB_IMG+sum(LSB_ANT)+sum(LSB_ANT_IMG)+sum(USB_ANT)+sum(USB_ANT_IMG);
%    sum_s=CSB+CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG);
   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
%  am_env=AM_signal;
 
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
 index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);
 
 
 
   y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,����30����λ��dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
     [~,id30]=max(FH(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
     
  [~,id9960]=max(FH(maxId+(9960-900)/freq_rev:maxId+(9960+900))/freq_rev);
   a9960=FH(maxId+(9960-900)/freq_rev-1+id9960);
                                  
                                     
                                     maxSuId_1=0;
                                      maxSuId_2=0;
                                      
                                      
                                    for iid=maxId+(9960-1000)/freq_rev:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_1=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    for iid=maxId+(9960+1000)/freq_rev:-1:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_2=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    ss=ifft(ft(maxId),NFFT);
                                   sss=real(ss);
                                   DC1=mean(sss);
                                   
                                   
                                    
                                     ss=ifft(ft(maxId:maxId+(30+210)/freq_rev),NFFT);
                                     am30_env=real(ss);
                                     am30_env=am30_env-mean(am30_env);
                                      
                                     AM30=2*max(am30_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(30,nfft,sss);
                                      
                                   
% F   =  [0:0.05:0.95]; 
% A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
% b  =  firls(20,F,A);
% 
%                                      Signal_Filter= filter(b,a,am30_env);
%                                      am30_env=Signal_Filter;
%                                      
%                                   
%                                    
%                                    am30_env=am30_env.*1000/5;
                                   
%                                      disp(['30HzAM: ',num2str(A_30/DC1)]);
                                     
%                                     sss=resample(yy,2^15,2^18);
%                                      fm30=lowp(abs(sss),40,90,1,30,2^15);%fft(sss,nfft);
%                                     
%                                      
%                                    am30_env=abs(fm30);
%                                      ft(maxId+(30+60)/freq_rev-1:maxId+(30+60)/freq_rev+1)=(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;

                                      ss=ifft(ft(maxId+(9960-480)/freq_rev:maxId+(9960+480)/freq_rev),NFFT);
                                   am9960_env=real(ss);
                                   
                                   AM9960=2*max(am9960_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(9960,nfft,sss);
%                                      disp(['9960HzAM: ',num2str(A_9960/DC1)]);
                                     
                                     
                                        diff9960=zeros(1,length(am9960_env));
                                     for i=1:length(am9960_env)-1
                                         diff9960(i)=(am9960_env(i+1)-am9960_env(i))/(1/fs);
                                     end
                                     
                                     sss=abs(hilbert(diff9960));


% sss=fmdemod(am9960_env,9960,rtlsdr_fs,480);    %diff(am9960_env);
                                    
%                                      sss=resample(rtl_fft,2^15,2^18);
%                                      fm30=lowp(sss,40,90,1,30,2^15);%fft(sss,nfft);
                                     fm30_env=real(sss)-mean(real(sss));
%                                      fm30_env(1)=fm30_env(2);
%                                      am30_env=resample(am30_env,2^15,nfft);
%                                      fm30_env=resample(fm30_env,2^15,nfft);
%                                      Wc=2*50/NFFT;                                          %��ֹƵ�� 50Hz
%                                      [b,a]=butter(4,Wc);
%                                      Signal_Filter=filter(b,a,fm30_env);
                                     
%                                      fm30_env=Signal_Filter;
                                     
                                        am30_env=resample(am30_env,1250,NFFT);
                                     fm30_env=resample(fm30_env,1250,NFFT);
                                     
                                       c_start=1;
                                   c_stop=length(am30_env)-c_start+1;  
                                     
                                   ww=blackman(c_stop-c_start+1);
%                                  ww=blackmanharris(c_stop-c_start+1);
                                     am30FFT=fft(am30_env(c_start:c_stop).*ww');
%                                     am30FFT=fft(am30_env(c_start:c_stop));
                                   am30FFTAMP=abs(am30FFT);
                                    [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                    id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30AM= angle(am30FFT(id30))*180/pi;
                                    
                                     fm30FFT=fft(fm30_env(c_start:c_stop).*ww');
%                                         fm30FFT=fft(fm30_env(c_start:c_stop));
                                   fm30FFTAMP=abs(fm30FFT);
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                   id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;
                                   

                                     R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
                                     [Rmax,Rloc]=max(R);
                                    Rloc=Rloc-(c_stop-c_start+1);
%                                     deg=Rloc*360*30/(rtlsdr_fs);

                                      deg=ph30FM-ph30AM;
                                      
                                  
            
                                    if deg<0
                                        deg=deg+360;
                                    end
                                    
                                    
 
                                    
                                    az_error=deg-ang;  %���㷽λ���
                                    
                                    az_error=az_error-180; 
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
                                    
                                  
                             
                                        
                                    
%                                     figure(2);
%                                     subplot(311);
%                                     plot(am30_env);
%                                    title('30Hz AM�ź�');
%                                    subplot(312);
%                                    plot(fm30_env);
%                                     title('30Hz FM�ź�');
%                                     subplot(313);
%                                     plot(R);
%                                     title(['XCORR ','Max Rloc==',num2str(Rloc)]);
 
                                      fmi=(maxSuId_2- maxSuId_1)*freq_rev/2/30;
                                     
                                      

                                          
                                          VOR_AZ=num2str(az_error);
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     if step_select==1
                                         step_unit=start_range+new_d;
                                     else
                                         step_unit=start_H+new_d;
                                     end
                                         results(sim_step,:)=[step_unit,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����
% xxlable="ʱ��";        %get(handles.edit_XLable,'String');
% yylable="����";         %get(handles.edit_YLable,'String');
%   figure(3);
  
%  
%    
%    subplot(3,2,3);
%     ftitle="�ռ�������ź�";
% plot(t,AM_signal);
% % xlabel(xxlable);
% ylabel(yylable);
% title(ftitle);
%     t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
   

 
 
 figure(3);
   subplot(3,2,6);
    ftitle="AZ: "+num2str(ang)+" | AM30="+VOR_30HzAM+" | AM9960="+ VOR_9960HzAM+" | AZ_E_r_r="+VOR_AZ;
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
 
 
 
 
 
%   
%   str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
%  t_p=0.07;
%  
%  text(t_p/2,1, str1,'Color','red','FontSize',8);
%  str2={num2str(FH(1)),  num2str(AM30_MOD,'%1.4f'),  num2str(AM9960_MOD,'%1.4f')};
%   text(t_p/2+33*t_p/200,1,str2,'Color','red','FontSize',8);
% 
% 
%    figure(10);
%     ftitle="Ƶ��ͼ";
%   freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
%   YY=fftshift(FH);
%   plot(freqaxis,YY);
% 
% xlabel("Ƶ��");
% ylabel("����");
% title(ftitle);
% grid on
%  

 
 
 
 
 %%
 
 
 
%  waitbar(0,hwait,'0%');
time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
   hwait.Name=['Left: ',num2str(fix(time_left)),'s ', num2str(time_over),'s/��',num2str(PerStr),'%'];
% pause(0.00005);
 
end    %����ѭ��������
 close(hwait);
 
%  toc
 
% LSB_ANT1=y.*LSB;
% figure(4);
% for k=1:sbs_total
%     plot(t,USB_ANT(k,:));
% hold on;
% end
% plot(t,b_sin);
% plot(t,b_cos);
% 
% 
%  t_p=T;
%  axis([-t_p t_p -2*CSB_A 2*CSB_A]);
% grid on;

 %%%%% results(sim_step,:)=[ang,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%����������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
     
 figure(11);
 cc=results(:,1);
 subplot(5,1,1);
  plot(cc,results(:,3));
     ftitle="VOR��λ��� Vs. ���� ";
xlabel("Radial(m)");
ylabel("AZ error(deg)");
title(ftitle);

 subplot(5,1,2);
  plot(cc,results(:,4));
     ftitle="30HzAM Vs. ���� ";
xlabel("Radial(m)");
ylabel("30HzAM(%)");
title(ftitle);

 subplot(5,1,3);
  plot(cc,results(:,5));
     ftitle="9960HzAM Vs. ���� ";
xlabel("Radial(m)");
ylabel("9960HzAM(%)");
title(ftitle);

 subplot(5,1,4);
  plot(cc,results(:,2));
     ftitle="RF LEVEL Vs. ����";
xlabel("Radial(m)");
ylabel("RF LEVEL��dBm)");
title(ftitle);


 subplot(5,1,5);
  plot(cc,results(:,6));
     ftitle="��Ƶָ��FMI Vs. ���� ";
xlabel("Radial(m)");
ylabel("FMI");
title(ftitle);
    
    
    
    
end     % ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

    
    
    
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
%   disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',ErrorInfo.message+"�к�:"+string(ErrorInfo.stack(1).line)+"��������"+ErrorInfo.stack(1).name);
end

%%%%%%%%%%%%%%%%%%%-------DVOR�������������------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ж�ĳֱ����ƽ��Ľ����Ƿ����ϰ���ƽ���ϣ�����ǣ���InObstale=true,����Ϊfalse.
function InObstacle=Inplane(p1,p2,Pttable,h)    %����boolֵ
    % p1��p2��ֱ�߶ε�����㡣
    % Pttable�ǹ����ϰ���ƽ����ĸ����㣬��һ������ƽ��
    % �㷨�����ĸ�������ȡ3������ƽ̨���̣�p1��p2����ֱ�߷��̣�������ý�������x,y,z.���ж�x,y,z�Ƿ��������ϣ�
    %ͨ�����㽻����������������֮�����ж��Ƿ�������ƽ���ڡ�
    pp=get(h.checkbox_Obstacle,'Value');  %�Ƿ�����ϰ������
if pp      %���ѡ����зɻ�λ�÷���
    
    [aa,bb]=size(Pttable);
    inplane=false;
    
    for kk=1:aa
                    fly_x=p1(1);
                    fly_y=p1(2);
                    fly_z=p1(3);
                    CSB_X=p2(1);
                    CSB_Y=p2(2);
                    CSB_Z=p2(3);
                    planes=Pttable{kk,1};  %ȡ��4���������ꡣ����CELL����
                 [aaa,bbb]=size(planes);
                 if aaa==0
                   continue;  %�����ǿ�����
                 else
                 A=planes(:,1); %������,���ﲻ����CELL����������
                 B=planes(:,2);
                 C=planes(:,3);
                 D=planes(:,4);
                 end
                    
                    
                        distan=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);
                        dirVector=[fly_x-CSB_X,fly_y-CSB_Y,fly_z-CSB_Z]/distan;

                       BA=B-A;
                       BC=B-C;
                       Nor_obs=cross(BA,BC);   %ƽ��������������cross ��˵õ�ƽ̨��������
                       

%                         
%  ֱ֪��L����m��m1��m2��m3�����ҷ�������ΪVL��v1��v2��v3����
%ƽ��P����n��n1��n2��n3�����ҷ��߷�������ΪVP��vp1��vp2��vp3�������ֱ����ƽ��Ľ���O�����꣨x��y��z����
% % ��ֱ�߷���д�ɲ���������ʽ�����У�
% % % x = m1+ v1 * t
% % y = m2+ v2 * t (1)
% % z = m3+ v3 * t
% % ��ƽ�淽��д�ɵ㷨ʽ������ʽ�����У�
% % vp1 * (x �C n1) + vp2 * (y �C n2) + vp3 * (z �C n3) = 0 (2)
%t = ((n1-m1)*vp1+(n2-m2)*vp2+(n3-m3)*vp3) / (vp1* v1+ vp2* v2+ vp3* v3);
%�����3��ʽ�з�ĸ(vp1* v1+ vp2* v2+ vp3* v3)Ϊ0�����ʾֱ����ƽ��ƽ�У���ֱ����ƽ��û�н��㡣
%����t��Ȼ��t����ʽ��1��������ý���O�����꣨x��y��z����
                        
             m1=CSB_X;
             m2=CSB_Y;
             m3=CSB_Z;
             v1=dirVector(1);
                v2=dirVector(2);   
                  v3=dirVector(3);
                  n1=B(1);
                  n2=B(2);
                  n3=B(3);
                  vp1=Nor_obs(1);
                   vp2=Nor_obs(2);
                    vp3=Nor_obs(3);
                   if (vp1* v1+ vp2* v2+ vp3* v3) ==0
                        
                       inplane=false;
                       
                   else
                       pp=((n1-m1)*vp1+(n2-m2)*vp2+(n3-m3)*vp3) / (vp1* v1+ vp2* v2+ vp3* v3);
                       x=m1+v1*pp;
                       y=m2+v2*pp;
                       z=m3+v3*pp;
%                        plot3(x,y,z,'o','Color',[1,0,0],'LineWidth',2); 
                       Pt=[x; y; z];
                       d_AB=norm(cross(B-A,Pt-A))/norm(B-A);   %���㽻�㵽�����ߵľ��룬���ǹ�ʽ���ò�˵õ�������������ٳ��Աߣ��͵õ��㵽ֱ�ߵľ���
                        d_BC=norm(cross(C-B,Pt-B))/norm(C-B); 
                         d_CD=norm(cross(D-C,Pt-C))/norm(D-C); 
                          d_DA=norm(cross(A-D,Pt-D))/norm(A-D); 
                          %���4�����������֮�͵������ε��������Pt���������ڣ�������
                          Diff_AB=bsxfun(@minus,A,B);
                          Dist_AB=sqrt(sum(Diff_AB.^2,1));
                          
                          Diff_BC=bsxfun(@minus,B,C);
                          Dist_BC=sqrt(sum(Diff_BC.^2,1));
                          
                          Diff_CD=bsxfun(@minus,C,D);
                          Dist_CD=sqrt(sum(Diff_CD.^2,1));
                          
                          Diff_DA=bsxfun(@minus,D,A);
                          Dist_DA=sqrt(sum(Diff_DA.^2,1));
                          
                          Spt=Dist_AB*d_AB+Dist_BC*d_BC+Dist_CD*d_CD+Dist_DA*d_DA;
                          d_CD2AB=norm(cross(B-A,C-A))/norm(B-A);
                          Splane=(Dist_CD+Dist_AB)*d_CD2AB;
                          
                          if abs(Spt-Splane)<0.000001
                              inplane=true;
 
 break;
                             
                          else
                              inplane=false;
                              
                          end
                              
                          
                          
                       
                   end
                   
  
%        if InObstacle   %�����ĳ���ϰ���ƽ�ϲ������ڵ����򷵻أ��� ·�źŲ������ۼӡ�
%             return;
%        end
    end
     InObstacle=inplane; 
else
   InObstacle=false;
    
end
             %%%%%%%%%%%------�ж��ɻ����ز��������Ƿ�����ཻ


function txt_freq_Callback(hObject, eventdata, handles)
% hObject    handle to txt_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_freq as text
%        str2double(get(hObject,'String')) returns contents of txt_freq as a double


% --- Executes during object creation, after setting all properties.
function txt_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% gain = handles.gain;                                % set gain as part of structure
% middle_gain = gain;                                 % set middle value to gain
freq_value = get(hObject, 'Value');                 % obtain value from slider
freq_value_str = num2str(freq_value);        % convert to string for output
% gain_output = gain + (freq_value-middle_gain);      % compute the new position of the gain value
% gain_output_num = str2num(num2str(gain_output,'%.1f')); % limit output to 1 decimal place
set(handles.txt_freq,'string',freq_value_str);  % set output to show user current gain

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% maximum_freq = get(hObject,'Max');  % obtain max freq value
% minimum_freq=get(hObject,'Min');    % obtain min freq value
% 
% middle_freq = (maximum_freq-minimum_freq)/2;                       % obtain middle value from max
% handles.middle_freq = middle_freq;                  % set gain as part of structure
% guidata(hObject, handles);                          % add it to the guidata for global use
% initial_freq = handles.initial_gain;                % use initial gain value from structure
% set(hObject,'Value', initial_freq);   
%     
    
    




function txt_Counterpoint_H_Callback(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_Counterpoint_H as text
%        str2double(get(hObject,'String')) returns contents of txt_Counterpoint_H as a double


% --- Executes during object creation, after setting all properties.

function txt_Counterpoint_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_Counterpoint_R_Callback(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_Counterpoint_R as text
%        str2double(get(hObject,'String')) returns contents of txt_Counterpoint_R as a double


% --- Executes during object creation, after setting all properties.
function txt_Counterpoint_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on txt_Counterpoint_H and none of its controls.
function txt_Counterpoint_H_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_H (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyChar = eventdata.Character; %get(hObject,'Stringr');
  if (~alldigits(keyChar)) && isprintable(keyChar)
        beep;
%         currentText =get(hObject,'String');
%          strrep(currentText,keyChar,'');
%         set(hObject,'String',currentText);
        set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
%     pause(1);
         currentText =get(hObject,'String');
%          disp(currentText);
                 newText=strrep(currentText,keyChar,'');
%          disp(['New:' newText]);
          set(hObject,'String',newText);
%          pause(1);
       
          if ~alldigits( newText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
       return; 
  else
          currentText =get(hObject,'String');
          if ~alldigits( currentText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
  end
    guidata(hObject, handles);   %����gui���� 
%     if ~isempty(currentText) &&  alldigits(currentText)
%         set(hMessageLabel, 'String','Please enter a 5-digit zip-code')
%         jEditbox.setBorder(javax.swing.border.LineBorder(java.awt.Color.red, 3, false));
%          beep;
       
%         eventdata.Character='';
        
%         set(hObject,'
%         set(hMessageLabel, 'String','Invalid zip: should only contain digits')
%         jEditbox.setBorder(javax.swing.border.LineBorder(java.awt.Color.red, 3, false));
%     elseif length(currentText) ~= 5
%         set(hMessageLabel, 'String','Invalid zip: should have exactly 5 digits')
%         jEditbox.setBorder(javax.swing.border.LineBorder(java.awt.Color.red, 3, false));
%     else
%         set(hMessageLabel, 'String','');
%         jEditbox.setBorder(javax.swing.border.LineBorder(java.awt.Color.gray, 1, false));
        
   
  
  
  
  function flag = alldigits(str)
    flag = ~isnan(str2double(str));

% Check whether a character is printable
function flag = isprintable(keyChar)
    keyVal = double(keyChar);
    flag = ~isempty(keyVal) && keyVal > 31 && keyVal < 128;
    
    
    function setSBants(handles,ants,low_high)   %low_high��һ��ָʾ1-25��26-48.
   sbants=ants;
   sb_h=str2double(get(handles.SBs_H,'String'));
   ring_r=str2double(get(handles.edit_AntennaArray_Dimension,'String'));
   
   angle=360/sbants;
   
        if low_high==1
%             dt1=[ants/2,7];
             dt1=cell(sbants/2,7);
                for i=1:sbants/2
                    for j=1:7
                        dt1{i,1}=i;
                        dt1{i,2}=(i-1)*angle;
                        dt1{i,3}=ring_r/2;
                        dt1{i,4}=sb_h;
                        dt1{i,5}=100;
                        dt1{i,6}=0.00; % %3.3f
                        dt1{i,7}=true;
                    end
                end
                
%                 dt1={'1','0.0','6.5','1.3','100','0','1'}; % ...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true;...
%                '1','0.0','6.5','1.3','100','0',true};
%            whos dt1;
            set(handles.uitable1,'Data',dt1);
%              set(handles.uitable2,'Data',dt1);
            
        else
            dt1=cell(sbants/2,7);
                for i=1:sbants/2
                    for j=1:7
                        dt1{i,1}=i+sbants/2;
                        dt1{i,2}=(i-1+sbants/2)*angle;
                        dt1{i,3}=ring_r/2;
                        dt1{i,4}=sb_h;
                        dt1{i,5}=100;
                        dt1{i,6}=0.00; % %3.3f
                        dt1{i,7}=true;
                    end
                end
                

            set(handles.uitable2,'Data',dt1);
%              set(handles.uitable2,'Data',dt1);
            
        end
    

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over txt_Counterpoint_H.
function txt_Counterpoint_H_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on txt_freq and none of its controls.
function txt_freq_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to txt_freq (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyChar = eventdata.Character; %get(hObject,'Stringr');
  if (~alldigits(keyChar)) && isprintable(keyChar)
        beep;
%         currentText =get(hObject,'String');
%          strrep(currentText,keyChar,'');
%         set(hObject,'String',currentText);
        set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
%     pause(1);
         currentText =get(hObject,'String');
%          disp(currentText);
                 newText=strrep(currentText,keyChar,'');
%          disp(['New:' newText]);
          set(hObject,'String',newText);
%          pause(1);
       
          if ~alldigits( newText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
        
  else
          currentText =get(hObject,'String');
          if ~alldigits( currentText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
  end
 drawnow;


% --- Executes on key press with focus on txt_Counterpoint_R and none of its controls.
function txt_Counterpoint_R_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_R (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyChar = eventdata.Character; %get(hObject,'Stringr');
  if (~alldigits(keyChar)) && isprintable(keyChar)
        beep;
%         currentText =get(hObject,'String');
%          strrep(currentText,keyChar,'');
%         set(hObject,'String',currentText);
        set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
%     pause(1);
         currentText =get(hObject,'String');
%          disp(currentText);
                 newText=strrep(currentText,keyChar,'');
%          disp(['New:' newText]);
          set(hObject,'String',newText);
%          pause(1);
       
          if ~alldigits( newText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
        
  else
          currentText =get(hObject,'String');
          if ~alldigits( currentText)
          set(hObject,'BackgroundColor',[1.0 1.0 0.0]);
          else
              set(hObject,'BackgroundColor',[1.0 1.0 1.0]);  
          end
  end
 drawnow;



function CSB_H_Callback(hObject, eventdata, handles)
% hObject    handle to CSB_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CSB_H as text
%        str2double(get(hObject,'String')) returns contents of CSB_H as a double


% --- Executes during object creation, after setting all properties.
function CSB_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CSB_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on CSB_H and none of its controls.
function CSB_H_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to CSB_H (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function wave_Length_Callback(hObject, eventdata, handles)
% hObject    handle to wave_Length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave_Length as text
%        str2double(get(hObject,'String')) returns contents of wave_Length as a double


% --- Executes during object creation, after setting all properties.
function wave_Length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_Length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SBs_H_Callback(hObject, eventdata, handles)
% hObject    handle to SBs_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SBs_H as text
%        str2double(get(hObject,'String')) returns contents of SBs_H as a double


% --- Executes during object creation, after setting all properties.
function SBs_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SBs_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_default.
function pushbutton_default_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format long;

try
    % Multi-line editboxes are contained within a scroll-panel
    freq = get(handles.txt_freq,'String');
    freq_value=str2double(freq);
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
   set(handles.wave_Length,'String',num2str(wave_L));
   set(handles.CSB_H,'String',num2str(half_wave_L));
    set(handles.SBs_H,'String',num2str(half_wave_L));
   set(handles.txt_Error,'Visible','Off');
   set(handles.txt_Error,'String','');
   
   
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    
    
    end
catch
    % probably a single-line editbox
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
end

% --- Executes on button press in button_init_AntSystem.
function button_init_AntSystem_Callback(hObject, eventdata, handles)
% hObject    handle to button_init_AntSystem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format long;
try
    % Multi-line editboxes are contained within a scroll-panel
    freq = get(handles.txt_freq,'String');
    fmindex=get(handles.edit_FMI,'String');
    
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
    
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    
   set(handles.wave_Length,'String',num2str(wave_L));
   set(handles.CSB_H,'String',num2str(half_wave_L));
    set(handles.SBs_H,'String',num2str(half_wave_L));
    set(handles.edit_AntennaArray_Dimension,'String',num2str(array_d));
   set(handles.txt_Error,'Visible','Off');
   set(handles.txt_Error,'String','');
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
    setSBants(handles,sbants,1);   %��ʼ���ߴ����ߵ����������LOW
    setSBants(handles,sbants,0);   %��ʼ���ߴ����ߵ����������HIGH
     
     set(handles.checkbox_ODD_low,'Value',1);
       set(handles.checkbox_ODD_high,'Value',1);
         set(handles.checkbox_EVEN_low,'Value',1);
           set(handles.checkbox_EVEN_high,'Value',1);
             set(handles.checkbox_USB_ODD,'Value',1);
               set(handles.checkbox_USB_EVEN,'Value',1);
                 set(handles.checkbox_LSB_ODD,'Value',1);
                   set(handles.checkbox_LSB_EVEN,'Value',1);
                      
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
  disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
end


% --- Executes on button press in pushbutton_AntennaPattern_V.
function pushbutton_AntennaPattern_V_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AntennaPattern_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%------�����ߴ�ֱ�淽����ͼ------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
       
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������








  freq = get(handles.txt_freq,'String');
  
   fly_step=get(handles.edit_step,'String');
    simulate_step=str2double(fly_step)*pi/180;  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
   
    freq_value=str2double(freq);
    if ~isnan(freq_value)
    lamda=300/freq_value;
     beta=2.*pi/lamda; 
 n=2; 
%  step=0.0001;
 
 t=0.0: simulate_step:pi/2;
 d1=str2double(get(handles.CSB_H,'String')); 
 d2=str2double(get(handles.txt_Counterpoint_H,'String')); 
 R1=str2double(get(handles.txt_Counterpoint_R,'String')); 
 F1=zeros(1,length(t));
 
 
 for k=1:length(t)
     if(t(k)<=atan(d1/R1))
         d=d1+d2;
         
         
 W(k)=2*beta*d*sin(t(k));
 
%  z1=((n/2).*W)-n/2*beta*d;
%  z2=((1/2).*W)-1/2*beta*d;
 F1(k)=sin(W(k));
 
 F1(k)=sqrt(1+reflect_rate_ground^2-2*reflect_rate_ground*cos(W(k)));
 
 
     else
         d=d1;
          W(k)=2*beta*d*sin(t(k));
          
%  z1=((n/2).*W)-n/2*beta*d;
%  z2=((1/2).*W)-1/2*beta*d;
 F1(k)=sin(W(k));
  F1(k)=sqrt(1+reflect_rate_couterpoint^2-2*reflect_rate_couterpoint*cos(W(k)));
     end
 end
 
 K1=abs(F1);
 figure(2);
 fMax=max(K1);
 fMin=min(K1);
 
 zMax=[];
 zMin=[];
 up_arrow=1;
 for i=1:length(t)-1     %���Ҽ�ֵ
       if abs(K1(i))>fMax/2 && up_arrow
          if abs(K1(i+1))<abs(K1(i))
          zMax=[zMax t(i)];%[z x(i)];
          up_arrow=0;
          end
       end
         if abs(K1(i))<fMax/2 && ~up_arrow
            if abs(K1(i+1))>abs(K1(i))
          zMin=[zMin t(i)];%[z x(i)];
          up_arrow=1;
            end
         end
 end

 figure(2);
 subplot(1,2,1);
 
 polar(t,K1);
 hold on
  text(0.5,-0.5,"���ֵ="+zMax*180/pi+"��");
  text(-0.5,-0.5,"��Сֵ="+zMin*180/pi+"��");
  disp("���ֵ="+zMax*180/pi+"��");
  disp("��Сֵ="+zMin*180/pi+"��");
%   syms x
%      f=sin(beta*d*sin(x*pi/180));
% s=abs(diff(f));        %һ�׵���
% s2=diff(f,2);          %���׵���
% 
% h=real(double(solve(s)));%һ�׵���Ϊ��ĵ���ܾ��Ǽ�ֵ�㣬ע���ǿ��ܣ�������������α�
%  for i=1:length(h) 
%  if subs(s2,x,h(i))<0 
%      text(cos(h(i)*pi/180),sin(h(i)*pi/180),"���ֵ="+double(subs(f,x,h(i)))+",�Ƕ�="+num2str(h(i))+"��");
%  disp(['������' num2str(h(i)) '��ȡ�ü���ֵ������ֵΪ' num2str(double(subs(f,x,h(i))))]);
%  elseif subs(s2,x,h(i))>0 
%  disp(['������' num2str(h(i)) '��ȡ�ü�Сֵ����СֵΪ' num2str(double(subs(f,x,h(i))))]) ;
%  text(cos(h(i)*pi/180),sin(h(i)*pi/180),"��Сֵ="+double(subs(f,x,h(i)))+",�Ƕ�="+num2str(h(i))+"��");
%  else 
%  disp(['������' num2str(h(i)) '�����׵���ҲΪ0�����ڸõ㴦���������м���ֵ����Сֵ���޼�ֵ']);%%%����������α� 
%  end 
%  end
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
%     
%   figure(3);
%   fplot(@(x)abs(sin(beta*d*sin(x*pi/180))),[0,180]);
%     s3=fzero(@(x)abs(sin(beta*d*sin(x*pi/180))),80);
%    s4=roots(f);
%     disp(s3);
%     disp(s4);
% % [xmin]=fminbnd(f,0,360),
% % [fxmin]=double(subs(f,x,xmin))
% % y=-f(x);
% % [xmax]=fminbnd(y,0,360)
% % [fxmax]=double(subs(f,x,xmax))
% % fplot(f,[0,360]);hold on
% % plot([xmin,xmax],subs(f,x,[xmin,xmax]),'ro','LineWidth',5);
% % text(xmin,double(subs(f,x,xmin)+0.2),'��Сֵ=',num2str(fxmin));
% % text(xmax,double(subs(f,x,xmax)+0.2),'����ֵ=',num2str(fxmax));
% 



carrier_V_diagram(hObject, eventdata, handles);




% x = @(u,v) exp(-abs(u)/10).*sin(5*abs(v));
% y = @(u,v) exp(-abs(u)/10).*cos(5*abs(v));
% z = @(u,v) u;
% 
% subplot(1,2,2);
% fs = fsurf(x,y,z);
% fs.URange = [-30 60];
% fs.FaceAlpha = .5;





% --- Executes on button press in radiobutton_Circle.
function radiobutton_Circle_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Circle
v=get(hObject,'Value');
if v==1
    set(handles.text_R,'Visible','On');
     set(handles.edit_R,'Visible','On');
      set(handles.text21,'Visible','On');
       set(handles.text22,'Visible','On');
        set(handles.text23,'Visible','Off');
         set(handles.text24,'Visible','On');
         set(handles.text24,'String','��');
          set(handles.text25,'Visible','On');
          set(handles.text25,'String','��');
          
           set(handles.edit_S_R,'Visible','On');
            set(handles.edit_S_R,'String','0.0');
           
            set(handles.text_S_R,'Visible','On');
            set(handles.text_S_R,'String','��ʼ��λ');
            
             set(handles.text_E_R,'Visible','On');
             set(handles.text_E_R,'String','������λ');
             
              set(handles.edit_E_R,'Visible','On');
              set(handles.edit_E_R,'String','360.0');
              
               set(handles.text_S_H,'Visible','On');
                set(handles.edit_S_H,'Visible','On');
                 set(handles.text_E_H,'Visible','Off');
                  set(handles.edit_E_H,'Visible','Off');
                   set(handles.text_radial_deg,'Visible','Off');
                    set(handles.edit_Radial,'Visible','Off');
                    set(handles.edit_step,'String','0.1');
                     set(handles.pushbutton_Phasing,'Visible','Off');
                    
end
                

% --- Executes on button press in radiobutton_Radial.
function radiobutton_Radial_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Radial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Radial
v=get(hObject,'Value');
if v==1
    set(handles.text_R,'Visible','Off');
     set(handles.edit_R,'Visible','Off');
      set(handles.text21,'Visible','Off');
       set(handles.text22,'Visible','On');
        set(handles.text23,'Visible','On');
         set(handles.text24,'Visible','On');
         set(handles.text24,'String','��');
          set(handles.text25,'Visible','On');
          set(handles.text25,'String','��');
          
           set(handles.edit_S_R,'Visible','On');
           set(handles.edit_S_R,'String','5000');
           
            set(handles.text_S_R,'Visible','On');
            set(handles.text_S_R,'String','��ʼ����');
            
             set(handles.text_E_R,'Visible','On');
             set(handles.text_E_R,'String','��������');
             
             set(handles.edit_E_R,'Visible','On');
             set(handles.edit_E_R,'String','20000');
             
               set(handles.text_S_H,'Visible','On');
                set(handles.edit_S_H,'Visible','On');
                 set(handles.text_E_H,'Visible','On');
                  set(handles.edit_E_H,'Visible','On');
                   set(handles.text_radial_deg,'Visible','On');
                    set(handles.edit_Radial,'Visible','On');
                    set(handles.edit_step,'String','10');
                    set(handles.pushbutton_Phasing,'Visible','On');
                    
end


function edit_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_R as text
%        str2double(get(hObject,'String')) returns contents of edit_R as a double


% --- Executes during object creation, after setting all properties.
function edit_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S_R as text
%        str2double(get(hObject,'String')) returns contents of edit_S_R as a double


% --- Executes during object creation, after setting all properties.
function edit_S_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_S_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_E_R_Callback(hObject, eventdata, handles)
% hObject    handle to edit_E_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_E_R as text
%        str2double(get(hObject,'String')) returns contents of edit_E_R as a double


% --- Executes during object creation, after setting all properties.
function edit_E_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_E_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S_H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S_H as text
%        str2double(get(hObject,'String')) returns contents of edit_S_H as a double


% --- Executes during object creation, after setting all properties.
function edit_S_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_S_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_E_H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_E_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_E_H as text
%        str2double(get(hObject,'String')) returns contents of edit_E_H as a double


% --- Executes during object creation, after setting all properties.
function edit_E_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_E_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Radial_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Radial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Radial as text
%        str2double(get(hObject,'String')) returns contents of edit_Radial as a double


% --- Executes during object creation, after setting all properties.
function edit_Radial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Radial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_BlendingFunction.
function popupmenu_BlendingFunction_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_BlendingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_BlendingFunction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_BlendingFunction


% --- Executes during object creation, after setting all properties.
function popupmenu_BlendingFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_BlendingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_Ground_ReflectFactor_Callback(hObject, eventdata, handles)
% hObject    handle to txt_Ground_ReflectFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_Ground_ReflectFactor as text
%        str2double(get(hObject,'String')) returns contents of txt_Ground_ReflectFactor as a double


% --- Executes during object creation, after setting all properties.
function txt_Ground_ReflectFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_Ground_ReflectFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_Counterpoint_ReflectFactor_Callback(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_ReflectFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_Counterpoint_ReflectFactor as text
%        str2double(get(hObject,'String')) returns contents of txt_Counterpoint_ReflectFactor as a double


% --- Executes during object creation, after setting all properties.
function txt_Counterpoint_ReflectFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_ReflectFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Lock.
function pushbutton_Lock_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Lock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
title=get(hObject,'String');
if title=="����"
    set(hObject,'String','�ɱ༭');
    set(handles.uitable1,'ColumnEditable',[false,false,false,false,false,false,false]);
     set(handles.uitable2,'ColumnEditable',[false,false,false,false,false,false,false]);
%     set(handles.uitable2,'ColumnEditable',[true,true,true,true,true,true,true]);
else
      set(hObject,'String','����');
    set(handles.uitable1,'ColumnEditable',[true,true,true,true,true,true,true]);
    set(handles.uitable2,'ColumnEditable',[true,true,true,true,true,true,true]);
end
    
    


% --- Executes on button press in checkbox_ODD_low.
function checkbox_ODD_low_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ODD_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ODD_low
v=get(hObject,'Value');
amp=get(handles.edit_A,'String');
ph=get(handles.edit_P,'String');
csb_h=get(handles.CSB_H,'String');
sb_h=get(handles.SBs_H,'String');
d=get(handles.edit_AntennaArray_Dimension,'String');

htable = handle(handles.uitable1);
[a,b]=size(htable.Data);
if v==1

for i=1:a
    if mod(i,2)~=0
     htable.Data{i,3}=str2double(d)/2;  
     htable.Data{i,4}=str2double(sb_h); 
    htable.Data{i,5}=str2double(amp); 
     htable.Data{i,6}=str2double(ph); 
    htable.Data{i,7}=true;
    
    end
end
else     %�����ѡ�񣬾Ͳ���Ҫ�ı�A��P
   for i=1:a
       if mod(i,2)~=0
           
%     htable.Data{i,5}=str2double(amp); 
%      htable.Data{i,6}=str2double(ph);        
    htable.Data{i,7}=false;
    
       end
   end
end
    


% --- Executes on button press in checkbox_ODD_high.
function checkbox_ODD_high_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ODD_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ODD_high
v=get(hObject,'Value');
amp=get(handles.edit_A,'String');
ph=get(handles.edit_P,'String');
csb_h=get(handles.CSB_H,'String');
sb_h=get(handles.SBs_H,'String');
d=get(handles.edit_AntennaArray_Dimension,'String');

htable = handle(handles.uitable2);
[a,b]=size(htable.Data);
if v==1

for i=1:a
    if mod(i,2)~=0
     htable.Data{i,3}=str2double(d)/2;  
     htable.Data{i,4}=str2double(sb_h); 
    htable.Data{i,5}=str2double(amp); 
     htable.Data{i,6}=str2double(ph); 
    htable.Data{i,7}=true;
    
    end
end
else     %�����ѡ�񣬾Ͳ���Ҫ�ı�A��P
   for i=1:a
       if mod(i,2)~=0
           
%     htable.Data{i,5}=str2double(amp); 
%      htable.Data{i,6}=str2double(ph);        
    htable.Data{i,7}=false;
    
       end
   end
end

% --- Executes on button press in checkbox_EVEN_low.
function checkbox_EVEN_low_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_EVEN_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_EVEN_low
v=get(hObject,'Value');
amp=get(handles.edit_A,'String');
ph=get(handles.edit_P,'String');
csb_h=get(handles.CSB_H,'String');
sb_h=get(handles.SBs_H,'String');
d=get(handles.edit_AntennaArray_Dimension,'String');

htable = handle(handles.uitable1);
[a,b]=size(htable.Data);
if v==1

for i=1:a
    if mod(i,2)==0
     htable.Data{i,3}=str2double(d)/2;  
     htable.Data{i,4}=str2double(sb_h); 
    htable.Data{i,5}=str2double(amp); 
     htable.Data{i,6}=str2double(ph); 
    htable.Data{i,7}=true;
    
    end
end
else     %�����ѡ�񣬾Ͳ���Ҫ�ı�A��P
   for i=1:a
       if mod(i,2)==0
           
%     htable.Data{i,5}=str2double(amp); 
%      htable.Data{i,6}=str2double(ph);        
    htable.Data{i,7}=false;
    
       end
   end
end

% --- Executes on button press in checkbox_EVEN_high.
function checkbox_EVEN_high_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_EVEN_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_EVEN_high
v=get(hObject,'Value');
amp=get(handles.edit_A,'String');
ph=get(handles.edit_P,'String');
csb_h=get(handles.CSB_H,'String');
sb_h=get(handles.SBs_H,'String');
d=get(handles.edit_AntennaArray_Dimension,'String');

htable = handle(handles.uitable2);
[a,b]=size(htable.Data);
if v==1

for i=1:a
    if mod(i,2)==0
     htable.Data{i,3}=str2double(d)/2;  
     htable.Data{i,4}=str2double(sb_h); 
    htable.Data{i,5}=str2double(amp); 
     htable.Data{i,6}=str2double(ph); 
    htable.Data{i,7}=true;
    
    end
end
else     %�����ѡ�񣬾Ͳ���Ҫ�ı�A��P
   for i=1:a
       if mod(i,2)==0
           
%     htable.Data{i,5}=str2double(amp); 
%      htable.Data{i,6}=str2double(ph);        
    htable.Data{i,7}=false;
    
       end
   end
end

% --- Executes on button press in checkbox_USB_ODD.
function checkbox_USB_ODD_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_USB_ODD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_USB_ODD


% --- Executes on button press in checkbox_USB_EVEN.
function checkbox_USB_EVEN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_USB_EVEN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_USB_EVEN


% --- Executes on button press in checkbox_LSB_ODD.
function checkbox_LSB_ODD_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LSB_ODD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_LSB_ODD


% --- Executes on button press in checkbox_LSB_EVEN.
function checkbox_LSB_EVEN_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LSB_EVEN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_LSB_EVEN



function edit_FMI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FMI as text
%        str2double(get(hObject,'String')) returns contents of edit_FMI as a double


% --- Executes during object creation, after setting all properties.
function edit_FMI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AntennaArray_Dimension_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AntennaArray_Dimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AntennaArray_Dimension as text
%        str2double(get(hObject,'String')) returns contents of edit_AntennaArray_Dimension as a double


% --- Executes during object creation, after setting all properties.
function edit_AntennaArray_Dimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AntennaArray_Dimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Load_Pattern.
function pushbutton_Load_Pattern_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load_Pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,filepath]=uigetfile('.txt','��ѡ��ߴ�����ģ��');
tt=["���ߺ�","�Ƕ�(��)","����(��)","�߶�(��)","����(%)","��λ(��)","��Ч"];
if filename~=0
str=[filepath,filename];
fid=fopen(str,'r','n');  % -tģʽ����ʾ���ı���ʽ�������Ƕ�����ģʽ���ж�д
filespecs=[repmat('%s',1,7) ,'%[^\n\r]'];
dataArray = textscan(fid, filespecs,  'Delimiter', '\t', 'TextType', 'string',  'ReturnOnError', false);

fclose(fid);
 celldata=[dataArray{1:end-1}];
 
 [a b]=size(celldata);
    for i=2:a
     for j=1:b
         if j<b
             if ~isempty(celldata{i,j}) 
             dt2{i-1,j}=str2double(celldata(i,j));
             else
                 dt2{i-1,j}="";
             end
         else
            temp=celldata{i,j};
            if temp=="1"
                dt2{i-1,j}=true;
            else
                dt2{i-1,j}=false;
            end
         end
     end
    end
 [a,b]=size(dt2);   
  if mod(a,2)==0 % a��ż�� 
set(handles.uitable1,'Data',dt2(1:a/2,1:b));
set(handles.uitable2,'Data',dt2(a/2+1:end,1:b));
else % a������ end
set(handles.uitable1,'Data',dt2);
  end
end

% dt2=[a,b];

% for i=1:a
%     for j=1:b
%          if j<b
%              if ~isempty(dt{i,j}) && ~isnan(dt{i,j})
%              dt2(i,j)=string(dt(i,j));
%              else
%                  dt2(i,j)="";
%              end
%          else
%             temp=dt{i,j};
%             if temp
%                 dt2(i,j)="1";
%             else
%                 dt2(i,j)="0";
%             end
%          end
%     end
% end

      
%  dt2=cell2mat(dt);

% save uidata.mat dt2


% dlmwrite(str,dt2,'delimiter','\t','-append');


% --- Executes on button press in pushbutton_Save_Pattern.
function pushbutton_Save_Pattern_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save_Pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dt=get(handles.uitable1,'Data');
dt1=get(handles.uitable2,'Data');
[a,b]=size(dt);
% dt2=[a,b];

for i=1:a
    for j=1:b
       
         if j<b
             if ~isempty(dt{i,j}) && ~isnan(dt{i,j})
             dt2(i,j)=string(dt(i,j));
            
             else
                 dt2(i,j)="";
              
             end
         else
            temp=dt{i,j};
            if temp
                dt2(i,j)="1";
            else
                dt2(i,j)="0";
            end
         end
    end
end
for i=1:a
    for j=1:b
       
         if j<b
             if ~isempty(dt1{i,j}) && ~isnan(dt1{i,j})
             dt3(i,j)=string(dt1(i,j));
            
             else
                 dt3(i,j)="";
              
             end
         else
            temp=dt1{i,j};
            if temp
                dt3(i,j)="1";
            else
                dt3(i,j)="0";
            end
         end
    end
end

      
%  dt2=cell2mat(dt);

% save uidata.mat dt2

[filename,filepath]=uiputfile('.txt','save file');
tt=["���ߺ�","�Ƕ�(��)","����(��)","�߶�(��)","����(%)","��λ(��)","��Ч"];
if filename~=0
str=[filepath,filename];
fid=fopen(str,'wt');  % -tģʽ����ʾ���ı���ʽ�������Ƕ�����ģʽ���ж�д
fprintf(fid,[repmat('%s\t',1,size(dt2,2)),'\n'],tt'); %
fprintf(fid,[repmat('%s\t',1,size(dt2,2)),'\n'],dt2'); %
fprintf(fid,[repmat('%s\t',1,size(dt2,2)),'\n'],dt3'); %
fclose(fid);
% dlmwrite(str,dt2,'delimiter','\t','-append');
end

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
fmi_value = get(hObject, 'Value');                 % obtain value from slider
fmi_value_str = num2str(fmi_value);        % convert to string for output
% gain_output = gain + (freq_value-middle_gain);      % compute the new position of the gain value
% gain_output_num = str2num(num2str(gain_output,'%.1f')); % limit output to 1 decimal place
set(handles.edit_FMI,'string',fmi_value_str);  % set output to show user current gain

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_FMI_Invert.
function pushbutton_FMI_Invert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FMI_Invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format long;
try
    % Multi-line editboxes are contained within a scroll-panel
    freq = get(handles.txt_freq,'String');
    array_d=get(handles.edit_AntennaArray_Dimension,'String');
    
    freq_value=str2double(freq);
    a_d=str2double(array_d);
    
    if ~isnan(freq_value)
    wave_L=300/freq_value;
   
    fmindex=pi*a_d/wave_L;
    set(handles.edit_FMI,'String',num2str(fmindex));
   set(handles.txt_Error,'Visible','Off');
  
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','������ֵ����');
    end
    
   
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
  disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','������ֵ����');
end



function edit_A_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A as text
%        str2double(get(hObject,'String')) returns contents of edit_A as a double


% --- Executes during object creation, after setting all properties.
function edit_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_P_Callback(hObject, eventdata, handles)
% hObject    handle to edit_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_P as text
%        str2double(get(hObject,'String')) returns contents of edit_P as a double


% --- Executes during object creation, after setting all properties.
function edit_P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CSB_Power_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CSB_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CSB_Power as text
%        str2double(get(hObject,'String')) returns contents of edit_CSB_Power as a double


% --- Executes during object creation, after setting all properties.
function edit_CSB_Power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CSB_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LSB_AMP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LSB_AMP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LSB_AMP as text
%        str2double(get(hObject,'String')) returns contents of edit_LSB_AMP as a double


% --- Executes during object creation, after setting all properties.
function edit_LSB_AMP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LSB_AMP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LSB_Phase_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LSB_Phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LSB_Phase as text
%        str2double(get(hObject,'String')) returns contents of edit_LSB_Phase as a double


% --- Executes during object creation, after setting all properties.
function edit_LSB_Phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LSB_Phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LSB_USB_Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LSB_USB_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LSB_USB_Ratio as text
%        str2double(get(hObject,'String')) returns contents of edit_LSB_USB_Ratio as a double


% --- Executes during object creation, after setting all properties.
function edit_LSB_USB_Ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LSB_USB_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AM30_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AM30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AM30 as text
%        str2double(get(hObject,'String')) returns contents of edit_AM30 as a double


% --- Executes during object creation, after setting all properties.
function edit_AM30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AM30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Noise.
function checkbox_Noise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Noise


% --- Executes on key press with focus on txt_Counterpoint_ReflectFactor and none of its controls.
function txt_Counterpoint_ReflectFactor_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to txt_Counterpoint_ReflectFactor (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_AntennaPattern_H.
function pushbutton_AntennaPattern_H_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AntennaPattern_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step as text
%        str2double(get(hObject,'String')) returns contents of edit_step as a double


% --- Executes during object creation, after setting all properties.
function edit_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Phasing.
function pushbutton_Phasing_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Phasing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Ѱ�������λ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
global Ptable;
try
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
    
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
    
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
  fc1= freq_value*1e6;  %�ز�Ƶ�ʣ���λHz
  fc=freq_value*1e6;
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
    f_Lsb1=fc1-9960;       %LSBƵ��
  f_Usb1=fc1+9960;  
  omega=30;            % �ز�����Ƶ��
  ph30=0;              % ��Ƶ��λ
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
 fs=5e6;   %����Ƶ��
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%  noisy=randn(1,N+1);
 
  t=0:1/fs:T-1/fs;    %����2��ʱ��Ƭ���ڵ��ź�;
  csb_power=str2double(get(handles.edit_CSB_Power,'String'));     % CSB����,��λW
  Lsb_phase=str2double(get(handles.edit_LSB_Phase,'String'));     %LSB��λ����λ�ȡ�
  Lsb_Amp=str2double(get(handles.edit_LSB_AMP,'String'))/100;     %LSB����
  Usb_Amp=str2double(get(handles.edit_LSB_USB_Ratio,'String'))/100;   % USB����
  AM30=str2double(get(handles.edit_AM30,'String'))/100;      % 30Hz AM ���ƶ�
  CSB_A=sqrt(2*csb_power*R);      
  
  csb_sideband=AM30*sin(2*pi*omega*t+ph30*pi/180);
% clr_sideband_150=n*sin(2*pi*omega2*t);
   
%    FFT_30=sin(2*pi*omega*t);
   
csb_mod=csb_sideband;
  
  
%%%%%%%%%%------------������������-------------------
blending_f=get(handles.popupmenu_BlendingFunction,'Value');
w = 2*T_on;     %1/1440*2=1/720  �ߴ���ʱ�䣬Ϊ��������
t_sb=fix(w/T*length(t)); %�ߴ��򿪵�ʱ��Ƭ�����������ݳ��ȣ�
if (mod(t_sb,2))~=0     %ȡż��
    t_sb=t_sb+1;
end

shift_t=t_sb/2;    %ƽ��50%�����
shift_sb=zeros(1,shift_t); %����һ��ƽ�������,�����桢ż��

switch blending_f
    case 1   % 'COS^2/SIN^2'
      b_cos=cos(2*pi*(1/(4*T_on))*t).*cos(2*pi*(1/(4*T_on))*t);  %��Ϻ�����������T_on�������� T_on�Ǳߴ���ʱ��
      b_sin=sin(2*pi*(1/(4*T_on))*t).*sin(2*pi*(1/(4*T_on))*t); 
    case 2    %'COS/SIN'
        b_cos=cos(2*pi*(1/(4*T_on))*t);
         b_sin=sin(2*pi*(1/(4*T_on))*t);
%     case 'CSC'
        
    case 3   % 'SQUARE'
        b_cos=1;
        b_sin=1;
end

% w = 2*T_on;     %1/1440;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
%   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
% CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
LSB_func=2*pi*f_Lsb*t;
USB_func=2*pi*f_Usb*t;
CSB_func=2*pi*fc*t;

 % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% % % % % % % % % % % % % %   antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
% % % % % % % % % % % % % %     couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
% % % % % % % % % % % % % %     counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
% % % % % % % % % % % % % %   reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
% % % % % % % % % % % % % %   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
% % % % % % % % % % % % % %    CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
% % % % % % % % % % % % % %        SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
% % % % % % % % % % % % % %      FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
% % % % % % % % % % % % % %      start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %      end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %        start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ����
% % % % % % % % % % % % % %        end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %��ʼ����
% % % % % % % % % % % % % %    FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
% % % % % % % % % % % % % %     simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�

%%%%%%%%%%%%%%%%%%%%%%Բ�ܷ���%%%%%%%%%%%%%
start_phase=0;
stop_phase=360;  %���෶Χ
hwait=waitbar(0,'processing...0%','name','�ߴ��ź������У�please wait>>>');

sbs_total=T/T_on;
steps=sbs_total/100;
CSB_X=0;
CSB_Y=0;
CSB_Z=CSB_H;
CSB_elev=atan(CSB_H/counterpoint_R)*180/pi;   %CSB������Ե��������ǣ���Ϊ�ж����Է��������Ǵ����Ϊ����ο�

LSB_ANT     =zeros(sbants,length(t));
LSB_ANT_IMG=zeros(sbants,length(t));
USB_ANT     =zeros(sbants,length(t));
USB_ANT_IMG  =zeros(sbants,length(t));

sb_low=get(handles.uitable1,'Data');
sb_high=get(handles.uitable2,'Data');
sb_all=cat(1,sb_low,sb_high);
[a,b]=size(sb_all);
sb_basic_data=zeros(sbants,6);   %%1:X,2:Y,3:Z,4:A,5:P,6.ON/OFF  �ߴ����߻����������ݾ���

sb_basic_data_image_1_counterpoint=zeros(sbants,6);    %�Է������Ϊ��׼�ľ���
sb_basic_data_image_2_ground=zeros(sbants,6);          %�Դ��Ϊ��׼�ľ���

if ~isempty(sb_all{1,1})
for b=1:sbants     %���ɱߴ���������  1�����ߺţ�2���Ƕȣ�3�����룬4���߶ȣ�5�����ȣ�6����λ��7������
    sb_z=sb_all{b,4};
    sb_r=sb_all{b,3};
    sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}/100;
    sb_p=sb_all{b,6};
    sb_onoff=sb_all{b,7};
    
    sb_x=sb_r*sin(-sb_ang*pi/180);
    sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4};
%     sb_r=sb_all{b,3};     %ֻ��Z��ı䣬����X��Y������һ���ġ�
%     sb_ang=sb_all{b,2};
     sb_a=sb_all{b,5}*reflect_rate_couterpoint/100;
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data_image_1_counterpoint(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4}- 2*couterpoint_H;
%     sb_r=sb_all{b,3};
%     sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}*reflect_rate_ground/100;   %���Ƿ���ϵ��
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
  sb_basic_data_image_2_ground(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
end
else
    set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',"���ȳ�ʼ���ߴ����������");
    return;
    
end   
%%

k=2*pi/wave_L;   %��λ���������ڲ�����λ����



simulation_range=start_phase:simulate_step:stop_phase;
step_c=length(simulation_range);

results=zeros(step_c,6);

  

for sim_step=1:step_c
    tic;
    sb_phase=start_phase+(sim_step-1)*simulate_step;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=start_H-couterpoint_H;
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
    fly_angle=atan(fly_z/FlySimulate_Circle)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(FlySimulate_Radial*pi/180);
    fly_y=fly_d*cos(FlySimulate_Radial*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
          
%             tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    
      rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
        
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                   reflect_rate=reflect_rate_couterpoint  ;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                  
                   ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                    
              end
    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+ delta_csb_p;

  loss_csb=32.45+20*log10(fc1)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc1)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
                    %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
   CSB= csb_Loss*CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
        CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*(cos(CSB_func)*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
   end
  
    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�д�����
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
                 
 
   
  
    
  %  ����
  %  receiver_SB �ӳ�����8696��
[ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb,f_Usb,sb_phase,handles,counterpoint_R,hwait);
  
%%
% %  for i=1:sbs_total   %sbants*30  �ܵķ���ߴ��ź�����
%      
% %     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
% %     y=SB_ANT(i,:);
% %     [a,b_cos]=size(y);
% % for k=1:b_cos
% %     if y(k)<0
% %         y(k)=0;
% %     end
% % end
% %   
%      % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
% %   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
%     
%     if mod(i,2)==0      %ż������
%         
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%          
%       sbnum_LSB=mod(i,sbants);
%       if sbnum_LSB==0
%           sbnum_LSB=sbants;
%       end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%                    
%                       %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                           LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                              LSB_ANT(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                     
%                     
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R  
%               
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                   reflect_rate=reflect_rate_couterpoint  ;
%                else
%                
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                   reflect_rate=reflect_rate_ground  ;
%                     
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                           LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                              LSB_ANT_IMG(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                    
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,������USB����������      
%                     
%                      sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                          USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                               USB_ANT(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R  
%                     
%              
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                           reflect_rate=reflect_rate_couterpoint  ;
%                else
%                 
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                           reflect_rate=reflect_rate_ground  ;
%                
%                 end     %����ѡ���IF���
%                 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*reflect_rate*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%               
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                          USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                               USB_ANT_IMG(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                  
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%       
%               
%                
%           
% 
%     
%     else                    %%��������
%         
%         
%         
%        b_odd= rectpuls(t-fix(i/2)*w,w);
%        
%         
%       sbnum_LSB=mod(i,sbants);
%          if sbnum_LSB==0
%           sbnum_LSB=sbants;
%          end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                        %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                           LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                              LSB_ANT(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R  
%               
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                   reflect_rate=reflect_rate_couterpoint  ;
%                else
%                  
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                   reflect_rate=reflect_rate_ground  ;
%                    
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                       %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                            LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                              LSB_ANT_IMG(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ���� ��USB���侵��     
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                  
%                 
%                
%             
% %                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                     xx=sb_x+tt*(fly_x-sb_x);
%                     yy=sb_y+tt*(fly_y-sb_y);
%                     rr=sqrt(xx^2+yy^2);
%                       
%                       
%                if rr<=counterpoint_R  
%                     
%              
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                           reflect_rate=reflect_rate_couterpoint  ;
%                else
%                  
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                           reflect_rate=reflect_rate_ground  ;
%                     
%                 
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=sb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                                 USB_ANT_IMG(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                  
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%        
% %          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
% %          USB_ANT(i,:)=USB.*b_odd.*b_cos;
%     
%       
% 
%     
%     
%     
%     end
%     
%     PerStr=fix(i/steps);
%     waitstr=['processing.....',num2str(PerStr),'%'];
%     waitbar(i/sbs_total,hwait,waitstr);
% %     pause(0.0005);
%  end
%  %%
%  
 
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
%%   
  %%   %%%%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
   a=length(EVEN_ANT);
  b=length(ODD_ANT);
  if a>b
      ODD_ANT=[ODD_ANT zeros(1,a-b)];
  else
      EVEN_ANT=[EVEN_ANT zeros(1,b-a)];
  end
   sum_s=CSB+CSB_IMG+EVEN_ANT+ODD_ANT; 
%    sum_s=CSB+CSB_IMG+sum(LSB_ANT)+sum(LSB_ANT_IMG)+sum(USB_ANT)+sum(USB_ANT_IMG);
%    sum_s=CSB+CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG);
   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�

   y_data_dbm = 10*log10((FH.^2)/50);  %���㹦��dBmֵ��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
     [~,id30]=max(FH(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
     
  [~,id9960]=max(FH(maxId+(9960-900)/freq_rev:maxId+(9960+900))/freq_rev);
   a9960=FH(maxId+(9960-900)/freq_rev-1+id9960);
                                  
                                     
                                     maxSuId_1=0;
                                      maxSuId_2=0;
                                      
                                      
                                    for iid=maxId+(9960-1000)/freq_rev:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_1=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    for iid=maxId+(9960+1000)/freq_rev:-1:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_2=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    ss=ifft(ft(maxId),NFFT);
                                   sss=real(ss);
                                   DC1=mean(sss);
                                   
                                   
                                    
                                     ss=ifft(ft(maxId:maxId+(30+210)/freq_rev),NFFT);
                                     am30_env=real(ss);
                                     am30_env=am30_env-mean(am30_env);
                                      
                                     AM30=2*max(am30_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(30,nfft,sss);
                                      
                                   

                                     ft(maxId+(30+60)/freq_rev-1:maxId+(30+60)/freq_rev+1)=(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;

                                      ss=ifft(ft(maxId+(9960-6000)/freq_rev:maxId+(9960+6000)/freq_rev),NFFT);
                                   am9960_env=real(ss);
                                   
                                   AM9960=2*max(am9960_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(9960,nfft,sss);

                                     
                                     
                                        diff9960=zeros(1,length(am9960_env));
                                     for i=1:length(am9960_env)-1
                                         diff9960(i)=(am9960_env(i+1)-am9960_env(i))/(1/fs);
                                     end
                                     
                                     sss=abs(hilbert(diff9960));



                                     fm30_env=real(sss)-mean(real(sss));

                                     
                                        am30_env=resample(am30_env,1250,NFFT);
                                     fm30_env=resample(fm30_env,1250,NFFT);
                                     
                                       c_start=1;
                                   c_stop=length(am30_env)-c_start+1;  
                                     
                                   ww=blackman(c_stop-c_start+1);
%                                  ww=blackmanharris(c_stop-c_start+1);
                                     am30FFT=fft(am30_env(c_start:c_stop).*ww');
%                                     am30FFT=fft(am30_env(c_start:c_stop));
                                   am30FFTAMP=abs(am30FFT);
                                    [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                    id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30AM= angle(am30FFT(id30))*180/pi;
                                    
                                     fm30FFT=fft(fm30_env(c_start:c_stop).*ww');
%                                         fm30FFT=fft(fm30_env(c_start:c_stop));
                                   fm30FFTAMP=abs(fm30FFT);
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                   id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;
                                   

                                     R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
                                     [Rmax,Rloc]=max(R);
                                    Rloc=Rloc-(c_stop-c_start+1);
%                                     deg=Rloc*360*30/(rtlsdr_fs);

                                      deg=ph30FM-ph30AM;
                                      
                                   
            
                                    if deg<0
                                        deg=deg+360;
                                    end
                                    
                                    az_error=deg-FlySimulate_Radial;  %���㷽λ���
                                    
                                    az_error=az_error-180; 
 
                                  
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
                                    
                                    
                                    
                                      fmi=(maxSuId_2- maxSuId_1)*freq_rev/2/30;
                                     
                                      

                                          
                                          VOR_AZ=num2str(az_error);
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     results(sim_step,:)=[sb_phase,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����
% xxlable="ʱ��";        %get(handles.edit_XLable,'String');
% yylable="����";         %get(handles.edit_YLable,'String');
%   figure(3);

time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
    hwait.Name=['ʣ��ʱ��:',num2str(time_left),'�롣---',num2str(PerStr),'%'];
pause(0.00005);
 



end    %����ѭ��������



 close(hwait);
 


 %%%%% results(sim_step,:)=[sb_phase,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%������(��λ������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
     
 figure(11);
 cc=results(:,1);
 subplot(4,1,1);
  plot(cc,results(:,3));
     ftitle="VOR��λ��� Vs. �ߴ���λ";
xlabel("SB phase");
ylabel("AZ error");
title(ftitle);

 subplot(4,1,2);
  plot(cc,results(:,4));
     ftitle="30HzAM Vs. �ߴ���λ";
xlabel("SB phase");
ylabel("30HzAM");
title(ftitle);

 subplot(4,1,3);
  plot(cc,results(:,5));
     ftitle="9960HzAM Vs. �ߴ���λ";
xlabel("SB phase");
ylabel("9960HzAM");
title(ftitle);

 subplot(4,1,4);
  plot(cc,results(:,2));
     ftitle="RF LEVEL Vs. �ߴ���λ";
xlabel("SB phase");
ylabel("RF LEVEL��dBm)");
title(ftitle);


% 
% plot(t,LSB_ANT(2,:));
% hold on;
% 
% plot(t,LSB_ANT(24,:));
% hold on;
% plot(t,LSB_ANT(25,:));
% hold on;
% plot(t,LSB_ANT(26,:));
% hold on;
% 
% plot(t,LSB_ANT(45,:));
% hold on;
% plot(t,LSB_ANT(46,:));
% hold on;
% plot(t,LSB_ANT(47,:));
% hold on;
% plot(t,LSB_ANT(48,:));
% hold on;
% 
% plot(t,LSB_ANT(49,:));
% hold on;
% 
% plot(t,USB_ANT(48,:));
% hold on;
% plot(t,USB_ANT(49,:));






%   
%   [yf,tf]=rcosine(Fd,Fs,'sqrt',R,Delay); 
% N=80000;
% x=randint(10000,1);        %��������ź�128����
% x=x*2-1;                 %�����ź����ֻ�������  ˫������
% xt=zeros(1,80000);         %����Ҫ��ֵ������
% xt(1:8:end)=x;           %ÿ��8�����һ��x���е�ֵ
% clear x;
% y=filter(yf,tf,xt);      %xtͨ�������˲�����y���˲��������a
% clear xt;
% fs=1e10;             %����Ƶ��
% n=0:N-1;          
% dt1=1/fs;             
% tc=n*dt1;
% f=n*fs/N;
% 
% fc=2.5e9;              %�ز�Ƶ��
% m=cos(2*pi*fc*tc);    %�ز�
% 
% z=0.795*y.*m;              %AM�����ź�   z���Ƿ����Ҫ������ź�  25dbm
% 
% nfft=32768;                            %fft����
% window=hanning(32768);                 %����131072����ĺ�����                               
% noverlap=16384;                         %�ֶλ���ĵ���
% figure(1)
% [Pz,f]=psd(r,nfft,fs,window,noverlap);   %�����źŵĹ����� 
% plot(f,10*log10(Pz))
% grid on
    
    
    
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
%   disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',ErrorInfo.message+"�к�:"+string(ErrorInfo.stack(1).line)+"��������"+ErrorInfo.stack(1).name);
end

%%%%%%%%%%%%%%%%%%%-------DVOR�������������------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in checkbox_D_H.
function checkbox_D_H_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_D_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_D_H
    function  carrier_V_diagram(hObject, eventdata, handles)
% hObject    handle to button_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%--------DVOR����������------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Ptable;
format long;
try
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
    
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step)*pi/180;  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
    
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
%     half_wave_L=wave_L/2;
%     array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
  fc= freq_value*1e6;  %�ز�Ƶ�ʣ���λHz
  fc1=freq_value*1e6;
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
  omega=30;            % �ز�����Ƶ��
  ph30=0;              % ��Ƶ��λ
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
 fs=4e6;   %����Ƶ��
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%  noisy=randn(1,N+1);
 
  t=0:1/fs:T-1/fs;    %����2��ʱ��Ƭ���ڵ��ź�;
  csb_power=str2double(get(handles.edit_CSB_Power,'String'));     % CSB����,��λW
%   Lsb_phase=str2double(get(handles.edit_LSB_Phase,'String'));     %LSB��λ����λ�ȡ�
%   Lsb_Amp=str2double(get(handles.edit_LSB_AMP,'String'))/100;     %LSB����
%   Usb_Amp=str2double(get(handles.edit_LSB_USB_Ratio,'String'))/100;   % USB����
  AM30=str2double(get(handles.edit_AM30,'String'))/100;      % 30Hz AM ���ƶ�
  CSB_A=sqrt(2*csb_power*R);      
  
  CSB_X=0;
  CSB_Y=0;
  CSB_Z= CSB_H;
  k=2*pi/wave_L;   %��λ����
  
  
  
  csb_sideband=AM30*sin(2*pi*omega*t+ph30*pi/180);
 
   
csb_mod=csb_sideband;
  
%   CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
%   LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
%   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
  



CSB_func=2*pi*fc*t;


start_z=0;     %couterpoint_H;
start_d=sqrt(start_range^2);
start_ang=atan(start_z/start_d);

start_x=start_d*sin(FlySimulate_Radial*pi/180);
start_y=start_d*cos(FlySimulate_Radial*pi/180);

stop_z=end_H-couterpoint_H;
stop_d=sqrt(end_range^2-stop_z^2);
stop_x=stop_d*sin(FlySimulate_Radial*pi/180);
stop_y=stop_d*cos(FlySimulate_Radial*pi/180);

distan=sqrt((stop_z-start_z)^2+(stop_x-start_x)^2+(stop_y-start_y)^2);
dirVector=[stop_x-start_x,stop_y-start_y,stop_z-start_z]/distan;

step_select=get(handles.checkbox_D_H,'Value');  %�����Ļ�׼���Ǹ߶Ȼ��Ǿ���


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ����  ��Բ  �ز�ǿ�ȷֲ�  ����   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_sel=get(handles.radiobutton_Radial,'Value');

if mode_sel==1
            simulation_range=start_ang:simulate_step:pi/2;    %  �����Բ  
step_c=length(simulation_range);

results=zeros(step_c,2);
hwait=waitbar(0,'processing...0%','name','please wait>>>');

for sim_step=1:step_c
    tic;
    ang=(sim_step-1)*simulate_step; %�������� 
    new_d=start_d*cos(ang);
    fly_x=new_d*sin(FlySimulate_Radial*pi/180);
    fly_y=new_d*cos(FlySimulate_Radial*pi/180);
  
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�� ���� rad
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=start_d*sin(ang);
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
%     fly_d=sqrt(newpoint(1)^2+newpoint(2)^2+newpoint(3)^2);   %  FlySimulate_Circle^2-fly_z^2);
%     fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
%     fly_x=fly_d*sin(ang*pi/180);
%     fly_y=fly_d*cos(ang*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%               tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    
      rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                  reflect_rate=reflect_rate_couterpoint  ;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                   ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
              end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    
    if mod(delta_csb_p,pi)<0.05
        disp("������㣺"+ang*180/pi+"   ���̲"+ delta_csb_p*180/pi+"��");
    end
    csb_img_p=ph_fc*pi/180+ delta_csb_p;

  loss_csb=32.45+20*log10(fc1)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc1)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
             %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
 CSB=csb_Loss*CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
 CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*(cos(CSB_func)*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
  end
       
  
   
  
    
  %%  


%  for i=1:sbs_total   %sbants*30  �ܵķ���ߴ��ź�����
%      
% %     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
% %     y=SB_ANT(i,:);
% %     [a,b_cos]=size(y);
% % for k=1:b_cos
% %     if y(k)<0
% %         y(k)=0;
% %     end
% % end
% %   
%      % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
% %   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
%     
%     if mod(i,2)==0      %ż������
%         
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%          
%       sbnum_LSB=mod(i,sbants);
%       if sbnum_LSB==0
%           sbnum_LSB=sbants;
%       end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  loss_sb=32.45+20*log10(f_Lsb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                     LSB_ANT(i,:)=LSB_f.*b_even.*b_sin;
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%           
%               if rr<=counterpoint_R
%                 
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                else
% %                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Lsb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                     LSB_ANT_IMG(i,:)=LSB_f.*b_even.*b_sin;
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,������USB����������      
%                     
%                      sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Usb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                USB_ANT(i,:)=USB_f.*b_even.*b_sin;
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                                 
%                           sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                else
% %                    if  fly_angle<sb_elev_min_usb
%                           sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
% %                        
% %                    end
%                end     %����ѡ���IF���
%                 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                USB_ANT_IMG(i,:)=USB_f.*b_even.*b_sin;
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%       
%               
%                
%           
% 
%     
%     else                    %%%%%%%%      ------��������
%           
%        b_odd= rectpuls(t-fix(i/2)*w,w);
%        
%         
%       sbnum_LSB=mod(i,sbants);
%          if sbnum_LSB==0
%           sbnum_LSB=sbants;
%          end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                     LSB_ANT(i,:)=LSB_f.*b_odd.*b_cos;
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%               else
%                  
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% % %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                     LSB_ANT_IMG(i,:)=LSB_f.*b_odd.*b_cos;
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ���� ��USB���侵��     
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                USB_ANT(i,:)=USB_f.*b_odd.*b_cos;
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                      
% %                     if sbnum_LSB<=sbants/2
% %                     sbnum_USB=sbnum_LSB+sbants/2;
% %                     else
% %                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
% %                     end
%                      tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                 
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%               else
%                   
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
% %                        
% %                    end
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                USB_ANT_IMG(i,:)=USB_f.*b_odd.*b_cos;
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%        
% %          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
% %          USB_ANT(i,:)=USB.*b_odd.*b_cos;
%     
%       
% 
%     
%     
%     
%     end
%     
%     PerStr=fix(i/steps);
%     waitstr=['processing.....',num2str(PerStr),'%'];
%     waitbar(i/sbs_total,hwait,waitstr);
% %     pause(0.0005);
%  end
 %%
 
 
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%           ��Ƶ��ƽ        %%%%%%%%%%%%
%%   
  
  
   sum_s=CSB+CSB_IMG ; %+sum(LSB_ANT)+sum(LSB_ANT_IMG)+sum(USB_ANT)+sum(USB_ANT_IMG);
%    sum_s=CSB+CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG);
   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
%  index_30=30/freq_rev+1;
%  index_9960=9960/freq_rev+1;
% %  AM30_MOD=FH(index_30)/FH(1);
% %  AM9960_MOD=FH(index_9960)/FH(1);
%  
 
 
%    y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,����30����λ��dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=FH(maxId);
%      RF_Level= num2str(rflevel);  
      
       results(sim_step,:)=[ang,rflevel];  %,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
     
   
 
 %%
 
 
 
%  waitbar(0,hwait,'0%');
time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
    hwait.Name=['ʣ��ʱ��:',num2str(time_left),'�롣---'];
    
%    PerStr=fix(i/steps);
    waitstr=['processing.....',num2str(PerStr),'%','--',num2str(ang*180/pi),'��'];
    waitbar( PerStr/100,hwait,waitstr);  
    
    
pause(0.00005);
 
end    %����ѭ��������
 close(hwait);
 


 %%%%% results(sim_step,:)=[ang,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%����������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
 F1=results(:,2);
  t=results(:,1);
 fMin=min(F1);
 


  K1=F1;  %+abs(fMin)+30;
   
 fMax=max(K1);
 
 zMax=[];
 zMin=[];
 up_arrow=1;
 for i=1:length(t)-1
       if abs(K1(i))>fMax/2 && up_arrow
          if K1(i+1)<K1(i)
          zMax=[zMax t(i)];%[z x(i)];
          up_arrow=0;
          end
       end
         if K1(i)<fMax/2 && ~up_arrow
            if K1(i+1)>K1(i)
          zMin=[zMin t(i)];%[z x(i)];
          up_arrow=1;
            end
         end
 end

 figure(2);
 subplot(1,2,2);
%  plot(t,results(:,2));
 fMax=fMax*1e4;
  polar(t,K1*1e4);
 
 hold on
  text(fMax/2,-fMax/2,"���ֵ="+zMax*180/pi+"��");
  text(-fMax/2,-fMax/2,"��Сֵ="+zMin*180/pi+"��");
  disp("���ֵ="+zMax*180/pi+"��");
  disp("��Сֵ="+zMin*180/pi+"��");
 
 
%  
%  figure(2);
%  cc=results(:,1);
%  subplot(2,1,2);
% 
%   polar(cc,results(:,2));
     ftitle="RF LEVEL(*1e4) Vs. ��ֱ����";
% xlabel("����");
% ylabel("RF LEVEL��dBm)");
title(ftitle);
   
 figure(4);
y_dBm=10*log10((results(:,2).^2)/50/2)+30;
  plot(t*180/pi,y_dBm);
    
         ftitle="RF LEVEL Vs. ��ֱ����";
 xlabel("����");
 ylabel("RF LEVEL��dBm)");
title(ftitle);
    
end     % ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

    
    
    
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
%   disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',ErrorInfo.message+"�к�:"+string(ErrorInfo.stack(1).line)+"��������"+ErrorInfo.stack(1).name);
end

%%%%%%%%%%%%%%%%%%%-------DVOR�������������--------


% --- Executes on button press in pushbutton_3D.
function pushbutton_3D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%--------��ά��ģ ����������------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;
 
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
    ph_ref=get(handles.edit_AZ_Align,'String');  %30Hz���࣬���ڷ�λУ��
    
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
       az_align=str2double(ph_ref);
       
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
  fc1= freq_value*1e6;  %�ز�Ƶ�ʣ���λHz
  fc=freq_value*1e6;
  
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
   f_Lsb1=fc1-9960;       %LSBƵ��
  f_Usb1=fc1+9960;       %USBƵ��
  omega=30;            % �ز�����Ƶ��
  ph30=az_align;              % ��Ƶ��λ
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
 fs=5e6;   %����Ƶ��
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%%%%%%%%%%%%%%%%%%%%%%Բ�ܷ���%%%%%%%%%%%%%
start_angle=start_range;
stop_angle=end_range;

sbs_total=T/T_on;
steps=sbs_total/100;
CSB_X=0;
CSB_Y=0;
CSB_Z=CSB_H;
CSB_elev=atan(CSB_H/counterpoint_R)*180/pi;   %CSB������Ե��������ǣ���Ϊ�ж����Է��������Ǵ����Ϊ����ο�

% LSB_ANT=zeros(sbants,length(t));
% LSB_ANT_IMG=zeros(sbants,length(t));
% USB_ANT=zeros(sbants,length(t));
% USB_ANT_IMG=zeros(sbants,length(t));

sb_low=get(handles.uitable1,'Data');
sb_high=get(handles.uitable2,'Data');
sb_all=cat(1,sb_low,sb_high);
[a,b]=size(sb_all);
sb_basic_data=zeros(sbants,6);   %%1:X,2:Y,3:Z,4:A,5:P,6.ON/OFF  �ߴ����߻����������ݾ���

sb_basic_data_image_1_counterpoint=zeros(sbants,6);    %�Է������Ϊ��׼�ľ���
sb_basic_data_image_2_ground=zeros(sbants,6);          %�Դ��Ϊ��׼�ľ���

if ~isempty(sb_all{1,1})
for b=1:sbants     %���ɱߴ���������  1�����ߺţ�2���Ƕȣ�3�����룬4���߶ȣ�5�����ȣ�6����λ��7������
    sb_z=sb_all{b,4};
    sb_r=sb_all{b,3};
    sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}/100;
    sb_p=sb_all{b,6};
    sb_onoff=sb_all{b,7};
    
    sb_x=sb_r*sin(-sb_ang*pi/180);
    sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4};
%     sb_r=sb_all{b,3};     %ֻ��Z��ı䣬����X��Y������һ���ġ�
%     sb_ang=sb_all{b,2};
     sb_a=sb_all{b,5}*reflect_rate_couterpoint/100;
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data_image_1_counterpoint(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4}-2*couterpoint_H;
%     sb_r=sb_all{b,3};
%     sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}*reflect_rate_ground/100;   %���Ƿ���ϵ��
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
  sb_basic_data_image_2_ground(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
end
else
    set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',"���ȳ�ʼ���ߴ����������");
    return;
    
end   
%%

DVOR3DModel(CSB_H,antenna_D,counterpoint_R,couterpoint_H,sb_basic_data,handles);  %�����ز����߷��ߣ�������ֱ���������뾶�������߶ȣ��ߴ����������





k=2*pi/wave_L;   %��λ���������ڲ�����λ����


%%%%%%%%%%%%%%%%%%%%%%%%%    Բ�ܷ��з���  ������������������������������

mode_sel=get(handles.radiobutton_Circle,'Value');
if mode_sel==1

simulation_range=start_angle:simulate_step:stop_angle;
step_c=length(simulation_range);

results=zeros(step_c,6);


for sim_step=1:step_c
    tic;
    ang=start_angle+(sim_step-1)*simulate_step;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=start_H-couterpoint_H;
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
    fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
%     tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    
      rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                  
                   ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
              end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+delta_csb_p;

    
  %%  


 for i=1:sbs_total   %sbants*30  �ܵķ���ߴ��ź�����
     

    
    if mod(i,2)==0      %ż������
        
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
         
      sbnum_LSB=mod(i,sbants);
      if sbnum_LSB==0
          sbnum_LSB=sbants;
      end
      sb_valid=sb_basic_data(sbnum_LSB,6);
      
      if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
                  sb_x=sb_basic_data(sbnum_LSB,1);
                  sb_y=sb_basic_data(sbnum_LSB,2);
                  sb_z=sb_basic_data(sbnum_LSB,3);
                  sb_a=sb_basic_data(sbnum_LSB,4);
                  sb_p=sb_basic_data(sbnum_LSB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
                    
                    sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
                    sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
                    
                   
 
                    
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
                    tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   %sb_elev_max_lsb
                sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
               else
%                    if  fly_angle<sb_elev_min_lsb
                        sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
                       sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
                      sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
                      sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
                      sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%              
                end
               
               
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
 
      else
           
      end     %����������ã�LSB��IF���
                        
                     
                    if sbnum_LSB<=sbants/2
                    sbnum_USB=sbnum_LSB+sbants/2;
                    else
                       sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
                    end
               
                
      sb_valid=sb_basic_data(sbnum_USB,6);
      
      if sb_valid   %����������ã��������Ӧ����,������USB����������      
                    
                     sb_x=sb_basic_data(sbnum_USB,1);
                  sb_y=sb_basic_data(sbnum_USB,2);
                  sb_z=sb_basic_data(sbnum_USB,3);
                  sb_a=sb_basic_data(sbnum_USB,4);
                  sb_p=sb_basic_data(sbnum_USB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
                    sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
                    sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                  
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
               tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
            
                    
                if rr<= counterpoint_R  %sb_elev_max_usb
                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
                         sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
                          sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
                          sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
                          sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
               else
%                    if  fly_angle<sb_elev_min_usb
                        sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
                          sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
                          sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
                          sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
                          sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                     
%                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
%                        
%                    end
                end     %����ѡ���IF���
                
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                   
      else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
      end    %�����������
      
              
               
          

    
    else                    %%��������
        
        
        
%        b_odd= rectpuls(t-fix(i/2)*w,w);
       
        
      sbnum_LSB=mod(i,sbants);
         if sbnum_LSB==0
          sbnum_LSB=sbants;
         end
      sb_valid=sb_basic_data(sbnum_LSB,6);
      
      if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
                  sb_x=sb_basic_data(sbnum_LSB,1);
                  sb_y=sb_basic_data(sbnum_LSB,2);
                  sb_z=sb_basic_data(sbnum_LSB,3);
                  sb_a=sb_basic_data(sbnum_LSB,4);
                  sb_p=sb_basic_data(sbnum_LSB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
                    
                    sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
                    sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
                    
                   
                    %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                   
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
                         tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   %sb_elev_max_lsb
                sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
               else
%                  
                        sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%      
                end
                      
%              
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
                    
 
                  
                    
      else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
          
      end     %����������ã�LSB��IF���
                        
                     
                    if sbnum_LSB<=sbants/2
                    sbnum_USB=sbnum_LSB+sbants/2;
                    else
                       sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
                    end
               
                
      sb_valid=sb_basic_data(sbnum_USB,6);
      
      if sb_valid   %����������ã��������Ӧ���� ��USB���侵��     
                    
                  sb_x=sb_basic_data(sbnum_USB,1);
                  sb_y=sb_basic_data(sbnum_USB,2);
                  sb_z=sb_basic_data(sbnum_USB,3);
                  sb_a=sb_basic_data(sbnum_USB,4);
                  sb_p=sb_basic_data(sbnum_USB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
                    sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
                    sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
 
                    
                
               
 
                       tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   
                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
                         sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
                          sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
                          sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
                          sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
               else      
%              
                        sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
                          sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
                          sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
                          sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
                          sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
                    
 
               end 
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                  
      else
% % %              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
% % %              USB_ANT_IMG(i,:)=zeros(1,length(t));
             
      end    %�����������
       
 
    
      

    
    
    
    end
    
 
 end
 
 
 
 
end    %����ѭ��������
 


end



  



function edit_AZ_Align_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AZ_Align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AZ_Align as text
%        str2double(get(hObject,'String')) returns contents of edit_AZ_Align as a double


% --- Executes during object creation, after setting all properties.
function edit_AZ_Align_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AZ_Align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_Obstacle.
function checkbox_Obstacle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Obstacle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Obstacle


% --- Executes on button press in pushbutton_Receiver.
function pushbutton_Receiver_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Receiver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%--------DVOR����������------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Ptable;

format long;
try
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
   ph_ref=get(handles.edit_AZ_Align,'String');  %30Hz���࣬���ڷ�λУ��
    
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
    az_align=str2double(ph_ref);
    
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
  fc1= freq_value*1e6;  %�ز�Ƶ�ʣ���λHz
  
  fc=freq_value*1e6;
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
  
   f_Lsb1=fc1-9960;       %LSBƵ��
  f_Usb1=fc1+9960; 
  
  omega=30;            % �ز�����Ƶ��
  ph30=az_align;              % ��Ƶ��λ
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
 fs=4e6;   %����Ƶ��,�������ԣ���4e6��6��8��10���Եõ��ԳƵ��ز�AM�ź� ����3��5e6��7��9���Գ�
% fs=500e3;
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%  noisy=randn(1,N+1);
 
  t=0:1/fs:T-1/fs;    %����0.1��ʱ��Ƭ���ڵ��ź�;
  csb_power=str2double(get(handles.edit_CSB_Power,'String'));     % CSB����,��λW
  Lsb_phase=str2double(get(handles.edit_LSB_Phase,'String'));     %LSB��λ����λ�ȡ�
  Lsb_Amp=str2double(get(handles.edit_LSB_AMP,'String'))/100;     %LSB����
  Usb_Amp=str2double(get(handles.edit_LSB_USB_Ratio,'String'))/100;   % USB����
  AM30=str2double(get(handles.edit_AM30,'String'))/100;      % 30Hz AM ���ƶ�
  CSB_A=sqrt(2*csb_power*R);      
  
  csb_sideband=AM30*sin(2*pi*omega*t+ph30*pi/180);
% clr_sideband_150=n*sin(2*pi*omega2*t);
   
%    FFT_30=sin(2*pi*omega*t);
   
csb_mod=csb_sideband;
  
  CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
  LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
  USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
  SB_signal=LSB+USB;
  sum_s=CSB+SB_signal;    %���ź�=�ز�+�ߴ�
  v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');   %��Ӹ�˹�������������źŵ�20dB
  
  else
   AM_signal=sum_s;
  end
  

  h=hilbert(AM_signal);   %����Hilbert�任������90��
%  yi=imag(h);       %Hilbert�任֮��õ��ĸ������鲿����Hilbert�任
%   xi=real(h);
 
  
  
  
 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 
 
 

 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
  
 
 ft=fft(am_env,NFFT);
 
 FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 
 
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
  
  index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 
 
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);

  
  
  
  ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����


xxlable="ʱ��";        %get(handles.edit_XLable,'String');
yylable="����";         %get(handles.edit_YLable,'String');
  figure(3);
  
  subplot(3,2,1);
plot(t,CSB);
% xlabel(xxlable);2
ylabel(yylable);
title(ftitle);
  t_p=T;
%  axis([0 t_p -CSB_A*(1+AM30) CSB_A*(1+AM30)]);
 
   subplot(3,2,2);
    ftitle="�ߴ�";
plot(t,SB_signal);
% xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -CSB_A CSB_A]);
   
   
   subplot(3,2,3);
    ftitle="�ռ�������ź�";
plot(t,AM_signal);
% xlabel(xxlable);
ylabel(yylable);
title(ftitle);
    t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);
   
    
   subplot(3,2,4);
    ftitle="AM���-----AM30="+num2str(AM30_MOD,'%1.4f')+" |  AM9960="+num2str(AM9960_MOD,'%1.4f');
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);
  
  str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
 t_p=0.07;
 
 text(t_p/2,1, str1,'Color','red','FontSize',8);
 str2={num2str(FH(1)),  num2str(FH(index_30),'%1.4f'),  num2str(FH(index_9960),'%1.4f')};
  text(t_p/2+33*t_p/200,1,str2,'Color','red','FontSize',8);
 
   subplot(3,2,5);
    ftitle="Ƶ��ͼ";
  freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
  YY=fftshift(FH);
  plot(freqaxis,YY);

xlabel("Ƶ��");
ylabel("����");
title(ftitle);
 
%  axis([-N N -0.2*CSB_A 2*CSB_A]);
grid on

%%%%%%%%%%------------������Ϻ��� ����-------------------
blending_f=get(handles.popupmenu_BlendingFunction,'Value');



switch blending_f
    case 1   % 'COS^2/SIN^2'
      b_cos=cos(2*pi*(1/(4*T_on))*t).*cos(2*pi*(1/(4*T_on))*t);  %��Ϻ�����������T_on�������� T_on�Ǳߴ���ʱ��
      b_sin=sin(2*pi*(1/(4*T_on))*t).*sin(2*pi*(1/(4*T_on))*t); 
    case 2    %'COS/SIN'
        b_cos=cos(2*pi*(1/(4*T_on))*t);
         b_sin=sin(2*pi*(1/(4*T_on))*t);
%     case 'CSC'
        
    case 3   % 'SQUARE'
        b_cos=1;
        b_sin=1;
end

  
w = 2*T_on;     %1/1440*2=1/720  �ߴ���ʱ�䣬Ϊ��������
t_sb=fix(w/T*length(t)); %�ߴ��򿪵�ʱ��Ƭ�����������ݳ��ȣ�
if (mod(t_sb,2))~=0     %ȡż��
    t_sb=t_sb+1;
end

shift_t=t_sb/2;    %ƽ��50%�����
shift_sb=zeros(1,shift_t); %����һ��ƽ�������,�����桢ż��

t_csb=t;

 
 
LSB_func=2*pi*f_Lsb*t;
USB_func=2*pi*f_Usb*t;
CSB_func=2*pi*fc*t;

 % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




% % % % % % % % % % % % % %   antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
% % % % % % % % % % % % % %     couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
% % % % % % % % % % % % % %     counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
% % % % % % % % % % % % % %   reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
% % % % % % % % % % % % % %   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
% % % % % % % % % % % % % %    CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
% % % % % % % % % % % % % %        SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
% % % % % % % % % % % % % %      FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
% % % % % % % % % % % % % %      start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %      end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %        start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ����
% % % % % % % % % % % % % %        end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %��ʼ����
% % % % % % % % % % % % % %    FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
% % % % % % % % % % % % % %     simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�

%%%%%%%%%%%%%%%%%%%%%%Բ�ܷ���%%%%%%%%%%%%%
start_angle=start_range;
stop_angle=end_range;
hwait=waitbar(0,'processing...0%','name','�ߴ��ź������У�please wait>>>');


sbs_total=T/T_on;
steps=sbs_total/100;
CSB_X=0;
CSB_Y=0;
CSB_Z=CSB_H;
CSB_elev=atan(CSB_H/counterpoint_R)*180/pi;   %CSB������Ե��������ǣ���Ϊ�ж����Է��������Ǵ����Ϊ����ο�

LSB_ANT=zeros(sbants,length(t));
LSB_ANT_IMG=zeros(sbants,length(t));
USB_ANT=zeros(sbants,length(t));
USB_ANT_IMG=zeros(sbants,length(t));

sb_low=get(handles.uitable1,'Data');
sb_high=get(handles.uitable2,'Data');
sb_all=cat(1,sb_low,sb_high);
[a,b]=size(sb_all);
sb_basic_data=zeros(sbants,6);   %%1:X,2:Y,3:Z,4:A,5:P,6.ON/OFF  �ߴ����߻����������ݾ���

sb_basic_data_image_1_counterpoint=zeros(sbants,6);    %�Է������Ϊ��׼�ľ���
sb_basic_data_image_2_ground=zeros(sbants,6);          %�Դ��Ϊ��׼�ľ���

if ~isempty(sb_all{1,1})
for b=1:sbants     %���ɱߴ���������  1�����ߺţ�2���Ƕȣ�3�����룬4���߶ȣ�5�����ȣ�6����λ��7������
    sb_z=sb_all{b,4};
    sb_r=sb_all{b,3};
    sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}/100;
    sb_p=sb_all{b,6};
    sb_onoff=sb_all{b,7};
    
    sb_x=sb_r*sin(-sb_ang*pi/180);
    sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4};
%     sb_r=sb_all{b,3};     %ֻ��Z��ı䣬����X��Y������һ���ġ�
%     sb_ang=sb_all{b,2};
     sb_a=sb_all{b,5}*reflect_rate_couterpoint/100;
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
    sb_basic_data_image_1_counterpoint(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
    sb_z=-sb_all{b,4}-2*couterpoint_H;
%     sb_r=sb_all{b,3};
%     sb_ang=sb_all{b,2};
    sb_a=sb_all{b,5}*reflect_rate_ground/100;   %���Ƿ���ϵ��
    sb_p=sb_all{b,6}+180;
%     sb_onoff=sb_all{b,7};
%     
%     sb_x=sb_r*sin(-sb_ang*pi/180);
%     sb_y=sb_r*cos(sb_ang*pi/180);
  sb_basic_data_image_2_ground(b,:)=[sb_x,sb_y,sb_z,sb_a,sb_p,sb_onoff];
    
end
else
    set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',"���ȳ�ʼ���ߴ����������");
    return;
    
end   
%%

k=2*pi/wave_L;   %��λ���������ڲ�����λ����


%%%%%%%%%%%%%%%%%%%%%%%%%    Բ�ܷ��з���  ������������������������������

mode_sel=get(handles.radiobutton_Circle,'Value');  %Բ�ܷ���ѡ��ť
if mode_sel==1       %ģʽѡ��=1����ʾԲ�ܷ���
simulation_range=start_angle:simulate_step:stop_angle;
step_c=length(simulation_range);

results=zeros(step_c,6);

step_c=1;   %����ѭ���������ջ��������ģʽ����������ͼ
for sim_step=1:step_c
    tic;
    ang=start_angle+(sim_step-1)*simulate_step;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=start_H-couterpoint_H;
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
    fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
%     tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);
    
      rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                  reflect_rate=reflect_rate_couterpoint  ;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                    ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
              end
   
     
    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+delta_csb_p;

  loss_csb=32.45+20*log10(fc1)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc1)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
  CSB= csb_Loss*CSB_A*(1+csb_mod).*cos(CSB_func+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
    CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*(cos(CSB_func)*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
   end
  
    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�д�����
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
    
  % receiver ģʽ
  %  receiver_SB �ӳ�����8696��
 [ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb,f_Usb,Lsb_phase,handles,counterpoint_R,hwait);
  

 %%
 
 
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
%%   
  %%   %%%%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
   a=length(EVEN_ANT);
  b=length(ODD_ANT);
  if a>b
      ODD_ANT=[ODD_ANT zeros(1,a-b)];
  else
      EVEN_ANT=[EVEN_ANT zeros(1,b-a)];
  end
   sum_s=CSB+CSB_IMG+EVEN_ANT+ODD_ANT;
%    sum_s=CSB+CSB_IMG+sum(LSB_ANT)+sum(LSB_ANT_IMG)+sum(USB_ANT)+sum(USB_ANT_IMG);
%     sum_s=CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG); %um(USB_ANT)+sum(USB_ANT_IMG);
   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 
 

 
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
 index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);
 ang30=angle(ft(index_30))*180/pi;
 
 
 
   y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,���ݷ��ȼ��㹦�ʣ���1/2 A^2/R,��λ��dBW,����30������dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
     [~,id30]=max(FH(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
     
  [~,id9960]=max(FH(maxId+(9960-900)/freq_rev:maxId+(9960+900))/freq_rev);
   a9960=FH(maxId+(9960-900)/freq_rev-1+id9960);
   AM9960=FH(maxId+(9960-900)/freq_rev+id9960-1)/FH(1)/(50/30);                               
                                     
                                     maxSuId_1=0;
                                      maxSuId_2=0;
                                      
                                      
                                    for iid=maxId+(9960-1000)/freq_rev:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_1=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    for iid=maxId+(9960+1000)/freq_rev:-1:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_2=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    ss=ifft(ft(maxId),NFFT);
                                   sss=real(ss);
                                   DC1=mean(sss);
                                   
                                   
                                    
                                     ss=ifft(ft(maxId:maxId+(30+210)/freq_rev),NFFT);
                                     am30_env=real(ss);
                                     am30_env=am30_env-mean(am30_env);
                                      
                                     AM30=2*max(am30_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(30,nfft,sss);
                                      
                                   
% F   =  [0:0.05:0.95]; 
% A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
% b  =  firls(20,F,A);
% 
%                                      Signal_Filter= filter(b,a,am30_env);
%                                      am30_env=Signal_Filter;
%                                      
%                                   
%                                    
%                                    am30_env=am30_env.*1000/5;
                                   
%                                      disp(['30HzAM: ',num2str(A_30/DC1)]);
                                     
%                                     sss=resample(yy,2^15,2^18);
%                                      fm30=lowp(abs(sss),40,90,1,30,2^15);%fft(sss,nfft);
%                                     
%                                      
%                                    am30_env=abs(fm30);

% 
%        ft(maxId+(30+60)/freq_rev-1:maxId+(30+60)/freq_rev+1)=(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;
% 
%                                       ss=ifft(ft(maxId+(9960-600)/freq_rev:maxId+(9960+600)/freq_rev),NFFT);
%                                    am9960_env=real(ss);







                                     ft(maxId+(30-20)/freq_rev:maxId+(30+100)/freq_rev)=0; %(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;

                                      ss=ifft(ft(maxId+(9960-6000)/freq_rev:maxId+(9960+6000)/freq_rev),NFFT);
%                                       ss=ifft(ft,NFFT);
                                   am9960_env=real(ss);
                             
                                   AM9960=2*max(am9960_env(5000/freq_rev:NFFT-5000/freq_rev))/DC1;%FindAmp(9960,nfft,sss);
%                                      disp(['9960HzAM: ',num2str(A_9960/DC1)]);
                                     
                                     
                                        diff9960=zeros(1,length(am9960_env));
                                     for i=1:length(am9960_env)-1
                                         diff9960(i)=(am9960_env(i+1)-am9960_env(i))/(1/fs);
                                     end
                                     
                                     sss=abs(hilbert(diff9960));


% sss=fmdemod(am9960_env,9960,rtlsdr_fs,480);    %diff(am9960_env);
                                    
%                                      sss=resample(rtl_fft,2^15,2^18);
%                                      fm30=lowp(sss,40,90,1,30,2^15);%fft(sss,nfft);
                                     fm30_env=real(sss)-mean(real(sss));
%                                      figure(10);
%                                      plot(filtfilt(b,a,fm30_env));
%                                      fm30_env(1)=fm30_env(2);
%                                      am30_env=resample(am30_env,2^15,nfft);
%                                      fm30_env=resample(fm30_env,2^15,nfft);
%                                      Wc=2*50/new_nfft;                                          %��ֹƵ�� 50Hz
%                                      [b,a]=butter(4,Wc);
%                                      Signal_Filter=filter(b,a,fm30_env);
%                                      
%                                      fm30_env=Signal_Filter;
                                     
                                        am30_env=resample(am30_env,1250,NFFT);
                                     fm30_env=resample(fm30_env,1250,NFFT);
                                     
                                       c_start=1;
                                   c_stop=length(am30_env)-c_start+1;  
                                     
                                   ww=blackman(c_stop-c_start+1);
%                                  ww=blackmanharris(c_stop-c_start+1);
                                     am30FFT=fft(am30_env(c_start:c_stop).*ww');
%                                     am30FFT=fft(am30_env(c_start:c_stop));
                                   am30FFTAMP=abs(am30FFT);
                                    [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                    id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30AM= angle(am30FFT(id30))*180/pi;
                                    
                                     fm30FFT=fft(fm30_env(c_start:c_stop).*ww');    %%%���ò�������������
%                                         fm30FFT=fft(fm30_env(c_start:c_stop));
                                   fm30FFTAMP=abs(fm30FFT);
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                   id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;
                                   

%                                      R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
%                                      [Rmax,Rloc]=max(R);
%                                     Rloc=Rloc-(c_stop-c_start+1);
%                                     deg=Rloc*360*30/1250;

                                      deg=ph30FM-ph30AM;
                                      
                                   
            
                                    if deg<0
                                        deg=deg+360;
                                    end
%                                     if deg<-180
%                                         deg=deg+360;
%                                     end
%                                     
%                                     ang_r=ang;
%                                     if ang>180
%                                         ang_r=ang-360;
%                                     end
                                    
                                    az_error=deg-ang;  %���㷽λ���
                                    
                                    az_error=az_error-180; 
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
%                                         if abs(az_error)>=180
%                                             az_error=az_error-360;
%                                         end
                                    
%                                     figure(2);
%                                     subplot(311);
%                                     plot(am30_env);
%                                    title('30Hz AM�ź�');
%                                    subplot(312);
%                                    plot(fm30_env);
%                                     title('30Hz FM�ź�');
%                                     subplot(313);
%                                     plot(R);
%                                     title(['XCORR ','Max Rloc==',num2str(Rloc)]);
 
                                      fmi=(maxSuId_2- maxSuId_1)*freq_rev/2/30;
                                     
                                      

                                          
                                          VOR_AZ=num2str(az_error);
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     results(sim_step,:)=[ang,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
 
 figure(3);
   subplot(3,2,6);
     ftitle="AZ: "+num2str(ang)+" | AM30="+VOR_30HzAM+" | AM9960="+ VOR_9960HzAM+" | AZ_E_r_r="+VOR_AZ;
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);

   
    figure(10);
   subplot(3,1,1);
    ftitle="AM������-AM30MOD";
plot(t,am_env);

xlabel(xxlable);
ylabel(yylable);
title(ftitle);


subplot(3,1,2);
    ftitle="9960AM���";
plot(t,am9960_env);

xlabel(xxlable);
ylabel(yylable);
title(ftitle);

subplot(3,1,3);
    ftitle="9960AM������";
plot(t,abs(hilbert(am9960_env)));

xlabel(xxlable);
ylabel(yylable);
title(ftitle);
 
time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
    hwait.Name=['Left: ',num2str(time_left),'s��---',num2str(PerStr),'%'];
pause(0.00005);
 
end    %����ѭ��������
 close(hwait);
 
 

 %%%%% results(sim_step,:)=[ang,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%����������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
     
 figure(11);
 cc=results(:,1);
 subplot(4,1,1);
  plot(cc,results(:,3));
     ftitle="VOR��λ��� Vs. �Ƕ�";
xlabel("Circle");
ylabel("AZ error");
title(ftitle);

 subplot(4,1,2);
  plot(cc,results(:,4));
     ftitle="30HzAM Vs. �Ƕ�";
xlabel("Circle");
ylabel("30HzAM");
title(ftitle);

 subplot(4,1,3);
  plot(cc,results(:,5));
     ftitle="9960HzAM Vs. �Ƕ�";
xlabel("Circle");
ylabel("9960HzAM");
title(ftitle);

 subplot(4,1,4);
  plot(cc,results(:,2));
     ftitle="RF LEVEL Vs. �Ƕ�";
xlabel("Բ��");
ylabel("RF LEVEL��dBm)");
title(ftitle);

end



 % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




% % % % % % % % % % % % % %   antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
% % % % % % % % % % % % % %     couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
% % % % % % % % % % % % % %     counterpoint_R=str2double(    D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
% % % % % % % % % % % % % %   reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
% % % % % % % % % % % % % %   reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
% % % % % % % % % % % % % %    CSB_H=str2double(   csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
% % % % % % % % % % % % % %        SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
% % % % % % % % % % % % % %      FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
% % % % % % % % % % % % % %      start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %      end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %��ʼ�߶�
% % % % % % % % % % % % % %        start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ����
% % % % % % % % % % % % % %        end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %��ʼ����
% % % % % % % % % % % % % %    FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
% % % % % % % % % % % % % %     simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
%������ʼ�����ά����
start_z=start_H-couterpoint_H;
start_d=sqrt(start_range^2-start_z^2);
start_x=start_d*sin(FlySimulate_Radial*pi/180);
start_y=start_d*cos(FlySimulate_Radial*pi/180);

stop_z=end_H-couterpoint_H;
stop_d=sqrt(end_range^2-stop_z^2);
stop_x=stop_d*sin(FlySimulate_Radial*pi/180);
stop_y=stop_d*cos(FlySimulate_Radial*pi/180);

distan=sqrt((stop_z-start_z)^2+(stop_x-start_x)^2+(stop_y-start_y)^2);
dirVector=[stop_x-start_x,stop_y-start_y,stop_z-start_z]/distan;

step_select=get(handles.checkbox_D_H,'Value');  %�����Ļ�׼���Ǹ߶Ȼ��Ǿ���


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    �������   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_sel=get(handles.radiobutton_Radial,'Value');  %����ѡ��ť�����ѡ�У������1

if mode_sel==1   %ѡ�������ģʽ��
            simulation_range=0:simulate_step:distan;
step_c=length(simulation_range);

results=zeros(step_c,6);
step_c=1;% ���ջ�ģʽ����ʾ���㲨��

for sim_step=1:step_c
    tic;
    new_d=(sim_step-1)*simulate_step; %�������� 
    
    newpoint=[start_x,start_y,start_z]+new_d.*dirVector; %�µ�����
    ang=FlySimulate_Radial;
    
    %simulation_range  %����ѭ����ʼ,�Ƕȵ�λ�Ƕȡ�
   %%�����ز����ߴ��ľ��룬���̲��λ����ȱ仯�������ɿռ��Ų���˥�����ù�ʽ
   %%L=32.45+20Lg(MHz)+20Lg(D)-GT(dB)-GR(dB),  Lת��Ϊ%�ȣ���CSB_A��ˣ��õ�Զ�˵ķ��ȡ�
   %%
    fly_z=newpoint(3)-couterpoint_H;
    if fly_z<0
       set(handles.txt_Error,'Visible','On');
       set(handles.txt_Error,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(newpoint(1)^2+newpoint(2)^2+newpoint(3)^2);   %  FlySimulate_Circle^2-fly_z^2);
    fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
    
    d_csb=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);  %�ز��Ĳ���
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%               tt=(0-(-CSB_Z))/(fly_z-(-CSB_Z));
%     xx=CSB_X+tt*(fly_x-CSB_X);
%     yy=CSB_Y+tt*(fly_y-CSB_Y);
%     rr=sqrt(xx^2+yy^2);

  rr=fly_d-fly_d*fly_z/(fly_z-(-CSB_Z));
    
    
                  A_x=CSB_X;
                  A_y=CSB_Y;
                  A_z=CSB_Z;
                  ZZ1=CSB_Z;
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز������ѡ��%%%%%%%%%%%%%%%%%
%                if fly_angle>CSB_elev
%               
%                   csb_z=-CSB_Z; 
%                  
%                else
%                    csb_z=-(CSB_Z+couterpoint_H);
%                end
              if rr<=counterpoint_R
                  csb_z=-CSB_Z;
                  reflect_rate=reflect_rate_couterpoint;
              else
                  csb_z=-(CSB_Z+2*couterpoint_H);
                    ZZ2=csb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-csb_z)/(fly_z-csb_z);
                   P_xx=A_x+tt*(fly_x-A_x);
                  P_yy=A_y+tt*(fly_y-A_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
              end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ ��2210�С�

    
    d_csb_image=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-csb_z)^2);  %�ز�����Ĳ���
    delta_d_csb=d_csb_image - d_csb;
    delta_csb_p=delta_d_csb*k+pi;   %�����ز����ߵ���λ��
    csb_img_p=ph_fc*pi/180+ delta_csb_p;

  loss_csb=32.45+20*log10(fc1)+20*log10(d_csb)-0-2.16-120-60;
                 csb_Loss=power(10,-loss_csb/20);
                 
                 loss_csb_img=32.45+20*log10(fc1)+20*log10(d_csb_image)-0-2.16-120-60;
                 csb_img_Loss=power(10,-loss_csb_img/20);
                 
              %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
    p1=[fly_x,fly_y,fly_z];
    p2=[CSB_X,CSB_Y,CSB_Z];
    p3=[CSB_X,CSB_Y,csb_z];
    checkCSB_OBS=Inplane(p1,p2,Ptable,handles);
     checkCSB_OBS_IMG=Inplane(p1,p3,Ptable,handles);
     
   if checkCSB_OBS   %���ڵ�
       CSB=zeros(1,length(t));
   else
  CSB= csb_Loss*CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
   end
   
   if checkCSB_OBS_IMG  %���ڵ�
       CSB_IMG=zeros(1,length(t));
   else
     CSB_IMG=csb_img_Loss*CSB_A*reflect_rate*(1+csb_mod).*(cos(CSB_func)*cos(csb_img_p)-sin(CSB_func)*sin(csb_img_p));
   end
  
    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�д�����
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
                 
     %  ������ջ�ģʽ������
    [ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb,f_Usb,Lsb_phase,handles,counterpoint_R,hwait);
  
    
  
    
  %%  


%  for i=1:sbs_total   %sbants*30  �ܵķ���ߴ��ź�����
%      
% %     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
% %     y=SB_ANT(i,:);
% %     [a,b_cos]=size(y);
% % for k=1:b_cos
% %     if y(k)<0
% %         y(k)=0;
% %     end
% % end
% %   
%      % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
% %   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
%     
%     if mod(i,2)==0      %ż������
%         
%          b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%          
%       sbnum_LSB=mod(i,sbants);
%       if sbnum_LSB==0
%           sbnum_LSB=sbants;
%       end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k,��λ�ǻ���
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
% 
%                       %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                             LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                   
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%           
%               if rr<=counterpoint_R
%                 
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                          reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                          reflect_rate=reflect_rate_ground;
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                            LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT_IMG(i,:)=LSB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                    
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,������USB����������      
%                     
%                      sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                  
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%              
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                     tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                                 
%                           sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                                  reflect_rate=reflect_rate_couterpoint;
%                else
% %                    if  fly_angle<sb_elev_min_usb
%                           sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                                  reflect_rate=reflect_rate_ground;
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
% %                        
% %                    end
%                end     %����ѡ���IF���
%                 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*reflect_rate*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT_IMG(i,:)=USB_f.*b_even.*b_sin;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%               
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%       
%               
%                
%           
% 
%     
%     else                    %%%%%%%%      ------��������
%           
%        b_odd= rectpuls(t-fix(i/2)*w,w);
%        
%         
%       sbnum_LSB=mod(i,sbants);
%          if sbnum_LSB==0
%           sbnum_LSB=sbants;
%          end
%       sb_valid=sb_basic_data(sbnum_LSB,6);
%       
%       if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
%                   sb_x=sb_basic_data(sbnum_LSB,1);
%                   sb_y=sb_basic_data(sbnum_LSB,2);
%                   sb_z=sb_basic_data(sbnum_LSB,3);
%                   sb_a=sb_basic_data(sbnum_LSB,4);
%                   sb_p=sb_basic_data(sbnum_LSB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     
%                    
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%   %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               LSB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                 LSB_ANT(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  
%                     
%                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                        tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                          reflect_rate=reflect_rate_couterpoint;
%               else
%                  
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                          reflect_rate=reflect_rate_ground;
% %                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
% % %                    end   
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                     LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*reflect_rate*(cos(LSB_func)*cos(phase_error)-sin(LSB_func)*sin(phase_error));
%                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                              LSB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                               LSB_ANT_IMG(i,:)=LSB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                     
%       else
%            LSB_ANT(i,:)=zeros(1,length(t));
%             LSB_ANT_IMG(i,:)=zeros(1,length(t));
%           
%       end     %����������ã�LSB��IF���
%                         
%                      
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
%                
%                 
%       sb_valid=sb_basic_data(sbnum_USB,6);
%       
%       if sb_valid   %����������ã��������Ӧ���� ��USB���侵��     
%                     
%                   sb_x=sb_basic_data(sbnum_USB,1);
%                   sb_y=sb_basic_data(sbnum_USB,2);
%                   sb_z=sb_basic_data(sbnum_USB,3);
%                   sb_a=sb_basic_data(sbnum_USB,4);
%                   sb_p=sb_basic_data(sbnum_USB,5);
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                    USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                               USB_ANT(i,:)=zeros(1,length(t));
%                        else
%                                USB_ANT(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%              
%                 
%                
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
%                      
% %                     if sbnum_LSB<=sbants/2
% %                     sbnum_USB=sbnum_LSB+sbants/2;
% %                     else
% %                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
% %                     end
%                      tt=(0-(-sb_z))/(fly_z-(-sb_z));
%                         xx=sb_x+tt*(fly_x-sb_x);
%                          yy=sb_y+tt*(fly_y-sb_y);
%                          rr=sqrt(xx^2+yy^2);
%     
%          
%               if rr<=counterpoint_R
%                 
%                         sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                                  reflect_rate=reflect_rate_couterpoint;
%               else
%                   
%                         sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
%                           sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
%                                  reflect_rate=reflect_rate_ground;
% %                     
% %                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
% %                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
% %                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
% %                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
% %                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
% %                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
% %                        
% %                    end
%                end 
%                  d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
%                     loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
%                  sb_Loss=power(10,-loss_sb/20);
%                  
%                     delta_phase=(d_sbant-d_csb)*k;   %���̲������λ����k
% %                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
% %                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
% %                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
%                     %    ֱ���ź� 
%                     phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
%                       USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*reflect_rate*sb_a*(cos(USB_func)*cos(phase_error)-sin(USB_func)*sin(phase_error));
%                     
%                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
%                         p1=[fly_x,fly_y,fly_z];
%                         p2=[sb_x,sb_y,sb_z];
%                         checkSB_OBS=Inplane(p1,p2,Ptable,handles);
%                        if checkSB_OBS   %���ڵ�
%                                USB_ANT_IMG(i,:)=zeros(1,length(t));
%                        else
%                                  USB_ANT_IMG(i,:)=USB_f.*b_odd.*b_cos;
%                        end
%                     %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %�д�����
%                      %
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    
%                 
%       else
%              USB_ANT(i,:)=zeros(1,length(t));   %������߲����ã�����0
%              USB_ANT_IMG(i,:)=zeros(1,length(t));
%              
%       end    %�����������
%        
% %          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
% %          USB_ANT(i,:)=USB.*b_odd.*b_cos;
%     
%       
% 
%     
%     
%     
%     end
%     
%     PerStr=fix(i/steps);
%     waitstr=['processing.....',num2str(PerStr),'%'];
%     waitbar(i/sbs_total,hwait,waitstr);
% %     pause(0.0005);
%  end
 %%
 
 
 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
%%   
  %%   %%%%%%%%%%%%%%%%%%%%%%%%%%������ز�����30Hz,9960Hz���ƶȣ���λ%%%%%%%%%%%%%
   a=length(EVEN_ANT);
  b=length(ODD_ANT);
  if a>b
      ODD_ANT=[ODD_ANT zeros(1,a-b)];
  else
      EVEN_ANT=[EVEN_ANT zeros(1,b-a)];
  end
   sum_s=CSB+CSB_IMG+EVEN_ANT+ODD_ANT;
%    sum_s=CSB+CSB_IMG+sum(LSB_ANT)+sum(LSB_ANT_IMG)+sum(USB_ANT)+sum(USB_ANT_IMG);
%    sum_s=CSB+CSB_IMG+sum(USB_ANT)+sum(USB_ANT_IMG);
   v=get(handles.checkbox_Noise,'Value');
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
  AM_signal=awgn(sum_s,20,'measured');
  
  else
   AM_signal=sum_s;
  end
  
      h=hilbert(AM_signal);   %����Hilbert�任������90��

 am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
%  am_env=AM_signal;
 
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
 index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);
 
 
 
   y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,����30����λ��dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
     [~,id30]=max(FH(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
     
  [~,id9960]=max(FH(maxId+(9960-900)/freq_rev:maxId+(9960+900))/freq_rev);
   a9960=FH(maxId+(9960-900)/freq_rev-1+id9960);
                                  
                                     
                                     maxSuId_1=0;
                                      maxSuId_2=0;
                                      
                                      
                                    for iid=maxId+(9960-1000)/freq_rev:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_1=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    for iid=maxId+(9960+1000)/freq_rev:-1:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_2=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    ss=ifft(ft(maxId),NFFT);
                                   sss=real(ss);
                                   DC1=mean(sss);
                                   
                                   
                                    
                                     ss=ifft(ft(maxId:maxId+(30+210)/freq_rev),NFFT);
                                     am30_env=real(ss);
                                     am30_env=am30_env-mean(am30_env);
                                      
                                     AM30=2*max(am30_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(30,nfft,sss);
                                      
                                   
% F   =  [0:0.05:0.95]; 
% A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
% b  =  firls(20,F,A);
% 
%                                      Signal_Filter= filter(b,a,am30_env);
%                                      am30_env=Signal_Filter;
%                                      
%                                   
%                                    
%                                    am30_env=am30_env.*1000/5;
                                   
%                                      disp(['30HzAM: ',num2str(A_30/DC1)]);
                                     
%                                     sss=resample(yy,2^15,2^18);
%                                      fm30=lowp(abs(sss),40,90,1,30,2^15);%fft(sss,nfft);
%                                     
%                                      
%                                    am30_env=abs(fm30);
                                     ft(maxId+(30+60)/freq_rev-1:maxId+(30+60)/freq_rev+1)=(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;

                                      ss=ifft(ft(maxId+(9960-6000)/freq_rev:maxId+(9960+6000)/freq_rev),NFFT);
                                   am9960_env=real(ss);
                                   
                                   AM9960=2*max(am9960_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(9960,nfft,sss);
%                                      disp(['9960HzAM: ',num2str(A_9960/DC1)]);
                                     
                                     
                                        diff9960=zeros(1,length(am9960_env));
                                     for i=1:length(am9960_env)-1
                                         diff9960(i)=(am9960_env(i+1)-am9960_env(i))/(1/fs);
                                     end
                                     
                                     sss=abs(hilbert(diff9960));


% sss=fmdemod(am9960_env,9960,rtlsdr_fs,480);    %diff(am9960_env);
                                    
%                                      sss=resample(rtl_fft,2^15,2^18);
%                                      fm30=lowp(sss,40,90,1,30,2^15);%fft(sss,nfft);
                                     fm30_env=real(sss)-mean(real(sss));
%                                      fm30_env(1)=fm30_env(2);
%                                      am30_env=resample(am30_env,2^15,nfft);
%                                      fm30_env=resample(fm30_env,2^15,nfft);
%                                      Wc=2*50/new_nfft;                                          %��ֹƵ�� 50Hz
%                                      [b,a]=butter(4,Wc);
%                                      Signal_Filter=filter(b,a,fm30_env);
%                                      
%                                      fm30_env=Signal_Filter;
                                     
                                        am30_env=resample(am30_env,1250,NFFT);
                                     fm30_env=resample(fm30_env,1250,NFFT);
                                     
                                       c_start=1;
                                   c_stop=length(am30_env)-c_start+1;  
                                     
                                   ww=blackman(c_stop-c_start+1);
%                                  ww=blackmanharris(c_stop-c_start+1);
                                     am30FFT=fft(am30_env(c_start:c_stop).*ww');
%                                     am30FFT=fft(am30_env(c_start:c_stop));
                                   am30FFTAMP=abs(am30FFT);
                                    [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                    id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30AM= angle(am30FFT(id30))*180/pi;
                                    
                                     fm30FFT=fft(fm30_env(c_start:c_stop).*ww');
%                                         fm30FFT=fft(fm30_env(c_start:c_stop));
                                   fm30FFTAMP=abs(fm30FFT);
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                   id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;
                                   

                                     R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
                                     [Rmax,Rloc]=max(R);
                                    Rloc=Rloc-(c_stop-c_start+1);
%                                     deg=Rloc*360*30/(rtlsdr_fs);

                                      deg=ph30FM-ph30AM;
                                      
                                  
            
                                    if deg<0
                                        deg=deg+360;
                                    end
                                    
                                    
 
                                    
                                    az_error=deg-ang;  %���㷽λ���
                                    
                                    az_error=az_error-180; 
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
                                    
                                  
                             
                                        
                                    
%                                     figure(2);
%                                     subplot(311);
%                                     plot(am30_env);
%                                    title('30Hz AM�ź�');
%                                    subplot(312);
%                                    plot(fm30_env);
%                                     title('30Hz FM�ź�');
%                                     subplot(313);
%                                     plot(R);
%                                     title(['XCORR ','Max Rloc==',num2str(Rloc)]);
 
                                      fmi=(maxSuId_2- maxSuId_1)*freq_rev/2/30;
                                     
                                      

                                          
                                          VOR_AZ=num2str(az_error);
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     if step_select==1
                                         step_unit=start_range+new_d;
                                     else
                                         step_unit=start_H+new_d;
                                     end
                                         results(sim_step,:)=[step_unit,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
      figure(3);
   subplot(3,2,6);
    ftitle="AM������-AM30MOD="+VOR_30HzAM+"  AM9960MOD="+ VOR_9960HzAM;
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);




  
    figure(10);
   subplot(3,1,1);
    ftitle="AM������-AM30MOD";
plot(am_env);

xlabel(xxlable);
ylabel(yylable);
title(ftitle);


subplot(3,1,2);
    ftitle="9960AM���";
plot(am9960_env);

xlabel(xxlable);
ylabel(yylable);
title(ftitle);

subplot(3,1,3);
    ftitle="9960AM������";
plot(abs(hilbert(am9960_env)));

xlabel(xxlable);
ylabel(yylable);
title(ftitle);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����
% xxlable="ʱ��";        %get(handles.edit_XLable,'String');
% yylable="����";         %get(handles.edit_YLable,'String');
%   figure(3);
  
%  
%    
%    subplot(3,2,3);
%     ftitle="�ռ�������ź�";
% plot(t,AM_signal);
% % xlabel(xxlable);
% ylabel(yylable);
% title(ftitle);
%     t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
   
    
%    subplot(3,2,6);
%     ftitle="AM������-AM30MOD="+num2str(AM30_MOD,'%1.4f')+"  AM9960MOD="+num2str(AM9960_MOD,'%1.4f');
% plot(t,am_env);
% xlabel(xxlable);
% ylabel(yylable);
% title(ftitle);
%    t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
%   
%   str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
%  t_p=0.07;
%  
%  text(t_p/2,1, str1,'Color','red','FontSize',8);
%  str2={num2str(FH(1)),  num2str(AM30_MOD,'%1.4f'),  num2str(AM9960_MOD,'%1.4f')};
%   text(t_p/2+33*t_p/200,1,str2,'Color','red','FontSize',8);
% 
% 
%    figure(10);
%     ftitle="Ƶ��ͼ";
%   freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
%   YY=fftshift(FH);
%   plot(freqaxis,YY);
% 
% xlabel("Ƶ��");
% ylabel("����");
% title(ftitle);
% grid on
%  

 
 
 
 
 %%
 
 
 
%  waitbar(0,hwait,'0%');
time_over=toc;
time_left=(step_c-sim_step)*time_over;
 PerStr=fix(sim_step/step_c*100);
    hwait.Name=['ʣ��ʱ��:',num2str(time_left),'�롣---',num2str(PerStr),'%'];
pause(0.00005);
 
end    %����ѭ��������
 close(hwait);
 
%%
%toc
 
% LSB_ANT1=y.*LSB;
% figure(4);
% for k=1:sbs_total
%     plot(t,USB_ANT(k,:));
% hold on;
% end
% plot(t,b_sin);
% plot(t,b_cos);
% 
% 
%  t_p=T;
%  axis([-t_p t_p -2*CSB_A 2*CSB_A]);
% grid on;
%%
 %%%%% results(sim_step,:)=[ang,rflevel,deg,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
 %%%%%%%%%%%%%%%%%%%%%%%%%%����������Ƶ��ƽ��vor��λ��30HzAM��9960HzAM��FMI
 %%%%%%%%%%%%%%%%%%%%%%%%%%  1       2        3        4        5       6       
     
 figure(11);
 cc=results(:,1);
 subplot(4,1,1);
  plot(cc,results(:,3));
     ftitle="VOR��λ��� Vs.������";
xlabel("distance");
ylabel("AZ error");
title(ftitle);

 subplot(4,1,2);
  plot(cc,results(:,4));
     ftitle="30HzAM Vs. ����";
xlabel("distance");
ylabel("30HzAM");
title(ftitle);

 subplot(4,1,3);
  plot(cc,results(:,5));
     ftitle="9960HzAM Vs. ����";
xlabel("distance");
ylabel("9960HzAM");
title(ftitle);

 subplot(4,1,4);
  plot(cc,results(:,2));
     ftitle="RF LEVEL Vs. ����";
xlabel("����");
ylabel("RF LEVEL��dBm)");
title(ftitle);
    
    
    
    
end     % ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

    
    
    
catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
%   disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',ErrorInfo.message+"�к�:"+string(ErrorInfo.stack(1).line)+"��������"+ErrorInfo.stack(1).name);
end

%%%%%%%%%%%%%%%%%%%-------DVOR�������������------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

    function [ODD_ANT, EVEN_ANT]=receiver_SB(t,t_sb,w,k,sbs_total,fly_x,fly_y,fly_z,b_cos,b_sin,Lsb_Amp,Usb_Amp,CSB_A,d_csb,reflect_rate_couterpoint,reflect_rate_ground,sbants,sb_basic_data,sb_basic_data_image_1_counterpoint,sb_basic_data_image_2_ground,LSB_func,USB_func,f_Lsb1,f_Usb1,Lsb_phase,handles,counterpoint_R,hwait)
            %%%%%%%%%%%%%%%%�ߴ�����$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            % t:����ʱ��,һ����0.1�롣 ��4e6�Ĳ����ʣ���t��4e5��������
            % t_sb����t����ʱ�䣨0.1�룩�ڣ���Ϻ����Ĳ�����Ϊż����
            % w: ��Ϻ����Ĵ򿪿��
            % k: �ز�Ƶ�ʵ���λ������ ��λ�Ƕ�/�׻򻡶�/�ס� k=2*pi/lambda
            %sbs_total: ��tʱ���ڣ��ߴ������źŵ������� 48*3Ȧ���ϱߴ���144�����±ߴ�Ҳ��144����
            %fly_x,y,z: �ɻ���������
            %b_cos��b_sin:  ��Ϻ���
            %Lsb_Amp: �±ߴ�����
            %Usb_Amp: �ϱߴ����ȣ�ʵ����һ�����±ߴ���صı���ֵ��
            %CSB_A: �ز�����
            %d_csb:  �ز�����
            %reflect_rate_couterpoint���������� ��
            %reflect_rate_ground ������淴�� �� 
            %sbants�� �ߴ����߸�����48��50
            %sb_basic_data���ߴ��������ݣ��ߴ�������*6����x,y,z����ͷ���A,��λp, ���ñ��
            %sb_basic_data_image_1_counterpoint�� �ߴ���������ڷ������ģ��Ļ�������
            %sb_basic_data_image_2_ground�� �ߴ�������Դ�صģ��Ļ�������
            %LSB_func��USB_func:�ϡ��±ߴ��ı��ʽ
            %f_Lsb1��f_Usb1: �ϡ��±ߴ���Ƶ��
            %Lsb_phase�� �±ߴ���λ
            %handles�����ھ��
            %counterpoint_R�������뾶�������жϷ����
            %hwait�����������
            
 global Ptable;   % actived obstacle list 
  
 
 wave_L=2*pi/k;
 freq_value=300/wave_L;
  fc=freq_value*1e6;
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
  
   
  
  wave_LSB=3e8/f_Lsb;
  wave_USB=3e8/f_Usb;
  
  k_LSB=2*pi/wave_LSB; 
  k_USB=2*pi/wave_USB;

 
 
EVEN_LSB=zeros(1,t_sb);
EVEN_USB=zeros(1,t_sb);

 
ODD_LSB=zeros(1,t_sb);
ODD_USB=zeros(1,t_sb); 
steps=sbs_total/100;
 for i=1:sbs_total   %    sbs_total=T/T_on  �ܵķ���ߴ��źţ��ڹ涨��ʱ��T�ڣ����� �ıߴ�������ÿ48��50һ��ѭ��
     
%     SB_ANT(i,:)=square(2*pi*30*(t-(i-1)*T_on),T_on/T*100);
%     y=SB_ANT(i,:);
%     [a,b_cos]=size(y);
% for k=1:b_cos
%     if y(k)<0
%         y(k)=0;
%     end
% end
%   
     % LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
%   USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
    
  if mod(i,2)==0      %ż������
        
         b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w); %���߿�ͨʱ�䣺һ����������
         for pp=1:length(b_even)-1
             if pp==1 && b_even(pp)==1
                 even_start=1;
             else
             if b_even(pp)==0 && b_even(pp+1)==1    %�ҳ�ÿ�����߷������ֹʱ��
                 even_start=pp+1;
             end
             end
             
             if pp==length(b_even)-1 && b_even(pp+1)==1
                 even_stop=pp+1;
             else
                if b_even(pp)==1 && b_even(pp+1)==0
                 even_stop=pp;
                end
             end
             
         end
           
             
          if b_sin==1
         b_even_sb=b_even(even_start:even_stop);
          else
         b_even_sb=b_even(even_start:even_stop).*b_sin(even_start:even_stop);  %��Ϻ���
          end 
         
% b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w).*b_sin*Lsb_Amp*CSB_A; %���߿�ͨʱ�䣺һ����������  

         
      sbnum_LSB=mod(i,sbants);
      if sbnum_LSB==0
          sbnum_LSB=sbants;    %�ߴ������1��ʼ������48��MOD����Ϊ0����ȡ48
      end
      sb_valid=sb_basic_data(sbnum_LSB,6);
      
      if sb_valid   %����������ã��������Ӧ����  �ȼ���LSB��������������
                  sb_x=sb_basic_data(sbnum_LSB,1);
                  sb_y=sb_basic_data(sbnum_LSB,2);
                  sb_z=sb_basic_data(sbnum_LSB,3);
                  
                 
                  
                  sb_a=sb_basic_data(sbnum_LSB,4);
                  sb_p=sb_basic_data(sbnum_LSB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                 
                    delta_phase=(d_sbant-d_csb)*k_LSB;   %���̲������λ����k,��λ�ǻ���
                    
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
                    
                   
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                    LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*cos(LSB_func(even_start:even_stop)+phase_error);

                    
                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                       if checkSB_OBS   %���ڵ�
                        EVEN_LSB=zeros(1,even_stop-even_start+1);
                           
%                         
                          
                       else
                           EVEN_LSB=LSB_f.*b_even_sb;
                           
                             
%                          
                       end
                    %%%%%%%%      �ϰ��ﷴ����棬         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                    
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
                      
                   A_x=sb_x;
                  A_y=sb_y;
                  A_z=sb_z;
                  ZZ1=sb_z;
                  
                    tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   %sb_elev_max_lsb
                sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
                  
                  
                  reflect_rate=reflect_rate_couterpoint;
             
                  
              
                  
               else
%                    if  fly_angle<sb_elev_min_lsb
                        sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
                       sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
                      sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
                        
                      sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
                      sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
                      
                      
                        ZZ2=sb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                    P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                                       
                       tt=(P_zz-sb_z)/(fly_z-sb_z);
                   P_xx=sb_x+tt*(fly_x-sb_x);
                  P_yy=sb_y+tt*(fly_y-sb_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
%                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                    end   
                end
               
               
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_LSB;   %���̲������λ����k
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                    LSB_f=sb_Loss*Lsb_Amp*CSB_A*reflect_rate*sb_a*cos(LSB_func(even_start:even_stop)+phase_error);

                      %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                         if checkSB_OBS   %���ڵ�
                             
                             EVEN_LSB=EVEN_LSB+zeros(1,even_stop-even_start+1);
%                         
                          
                         else
                              EVEN_LSB=EVEN_LSB+LSB_f.*b_even_sb;
                             
%                          
                        end
                        
                       
                    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                
                    
      else  %��� ���߲����ã���������档
                            
          EVEN_LSB=zeros(1,even_stop-even_start+1);
        
          
      end     %����������ã�LSB��IF���
      
    
                     
       %%%%%%%%%  ������USB    ************************************************************************      
      
      if sbnum_LSB<=sbants/2
                    sbnum_USB=sbnum_LSB+sbants/2;
                    else
                       sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
       end
               
                
      sb_valid=sb_basic_data(sbnum_USB,6);
      
      if sb_valid   %����������ã��������Ӧ����,������USB����������      
                    
                  sb_x=sb_basic_data(sbnum_USB,1);
                  sb_y=sb_basic_data(sbnum_USB,2);
                  sb_z=sb_basic_data(sbnum_USB,3);
                  sb_a=sb_basic_data(sbnum_USB,4);
                  sb_p=sb_basic_data(sbnum_USB,5);
                 
                  
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_USB;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                   USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*cos(USB_func(even_start:even_stop)+phase_error);                  
                            
                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                      if checkSB_OBS   %���ڵ�
                               
                             EVEN_USB=zeros(1,even_stop-even_start+1);
%                         
                          
                         else
                              EVEN_USB=USB_f.*b_even_sb;
                             
%                          
                       end
                        
                        
                    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
                   
              
                
               
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
              
                A_x=sb_x;
                  A_y=sb_y;
                  A_z=sb_z;
                  ZZ1=sb_z;
                  
               tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
            
                    
                if rr<= counterpoint_R  %sb_elev_max_usb
                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
                         sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
                          sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
                          sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
                          sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
                          reflect_rate=reflect_rate_couterpoint;
                          
               else
%                    if  fly_angle<sb_elev_min_usb
                        sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
                          sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
                          sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
                          sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
                          sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
                          
                          
                        ZZ2=sb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                          P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                          
                       tt=(P_zz-sb_z)/(fly_z-sb_z);
                   P_xx=sb_x+tt*(fly_x-sb_x);
                  P_yy=sb_y+tt*(fly_y-sb_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                          
                          
                          
                         
%                     
%                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5); %�ߴ����ߴ������߻���ϵͳ��·��������ơ�
%                        
%                    end
                end     %����ѡ���IF���
                
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                    delta_phase=(d_sbant-d_csb)*k_USB;   %���̲������λ����k
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %Ԥ�ñߴ���λ+������·��λ+���̲��λΪ����
                   USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A* reflect_rate*sb_a*cos(USB_func(even_start:even_stop)+phase_error); 
                   
                             
                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                       
                       if checkSB_OBS   %���ڵ�
                                
                             EVEN_USB=EVEN_USB+zeros(1,even_stop-even_start+1);
%                         
                          
                         else
                              EVEN_USB=EVEN_USB+USB_f.*b_even_sb;
                             
%                          
                        end 
                        
                    %%%%%%%%      ������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
            
      else   %��� ���߲����ã�������
             EVEN_USB=zeros(1,even_stop-even_start+1);
          
          
           
             
      end    %�����������
      
           if i==2
               EVEN_ANT=EVEN_LSB+EVEN_USB;
           else
               EVEN_ANT=[EVEN_ANT EVEN_LSB+EVEN_USB];
           end    
          

 %%   ODD�������Σ�
    
    else                    %%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
%        b_odd= rectpuls(t,w);
         b_odd= rectpuls(t-fix(i/2)*w,w);
         
         
          for pp=1:length(b_odd)-1
              if pp==1 && b_odd(pp)==1
                  odd_start=1;
              else
                  
             if b_odd(pp)==0 && b_odd(pp+1)==1    %�ҳ�ÿ�����߷������ֹʱ��
                 odd_start=pp+1;
             end
              end
              
              if pp==length(b_odd)-1 && b_odd(pp+1)==1
                  odd_stop=pp+1;
              else
               if b_odd(pp)==1 && b_odd(pp+1)==0  
                 odd_stop=pp;
               end
              end
             
          end
          if b_cos==1
            b_odd_sb=b_odd(odd_start:odd_stop);
          else
          b_odd_sb=b_odd(odd_start:odd_stop).*b_cos(odd_start:odd_stop);  %��Ϻ���
          end
         
         
      sbnum_LSB=mod(i,sbants);
         if sbnum_LSB==0
          sbnum_LSB=sbants;
         end
      sb_valid=sb_basic_data(sbnum_LSB,6);
      
      if sb_valid   %����������ã��������Ӧ����,����LSB���侵��
                  sb_x=sb_basic_data(sbnum_LSB,1);
                  sb_y=sb_basic_data(sbnum_LSB,2);
                  sb_z=sb_basic_data(sbnum_LSB,3);
                  sb_a=sb_basic_data(sbnum_LSB,4);
                  sb_p=sb_basic_data(sbnum_LSB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_LSB;   %���̲������λ����k
                    
%                     sb_elev_max_lsb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_lsb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
                    
                   
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                    LSB_f=sb_Loss*Lsb_Amp*CSB_A*sb_a*cos(LSB_func(odd_start:odd_stop)+phase_error);

                              
                     %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                        
                        if checkSB_OBS   %���ڵ�
                             
                         ODD_LSB=zeros(1,odd_stop-odd_start+1);
                           
%                         
                          
                       else
                           ODD_LSB=LSB_f.*b_odd_sb;
                             
%                          
                        end 
                        
                    %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                  
                    
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
                        A_x=sb_x;
                  A_y=sb_y;
                  A_z=sb_z;
                  ZZ1=sb_z;
                  
                      
                         tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   %sb_elev_max_lsb
                sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
                   reflect_rate=reflect_rate_couterpoint;
               else
%                    if  fly_angle<sb_elev_min_lsb
                        sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
                  sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
                  sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
                  sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
                  sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
                  
                  
                  
                        ZZ2=sb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                         P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                         
                       tt=(P_zz-sb_z)/(fly_z-sb_z);
                   P_xx=sb_x+tt*(fly_x-sb_x);
                  P_yy=sb_y+tt*(fly_y-sb_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                  
                  
                  
                  
%                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                    end   
                end
                      
%                if fly_angle>sb_elev_max_lsb
%                 sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                else
%                    if  fly_angle<sb_elev_min_lsb
%                         sb_x=sb_basic_data_image_2_ground(sbnum_LSB,1);
%                   sb_y=sb_basic_data_image_2_ground(sbnum_LSB,2);
%                   sb_z=sb_basic_data_image_2_ground(sbnum_LSB,3);
%                   sb_a=sb_basic_data_image_2_ground(sbnum_LSB,4);
%                   sb_p=sb_basic_data_image_2_ground(sbnum_LSB,5);
%                     else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_LSB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_LSB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_LSB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_LSB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_LSB,5);
%                    end   
%                end 
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Lsb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_LSB;   %���̲������λ����k
                    
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                    LSB_f=sb_Loss*Lsb_Amp*reflect_rate*CSB_A*sb_a*cos(LSB_func(odd_start:odd_stop)+phase_error);
  %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                        
                     if checkSB_OBS   %���ڵ�
                         
                               ODD_LSB=ODD_LSB+zeros(1,odd_stop-odd_start+1);
%                         
                          
                         else
                              ODD_LSB=ODD_LSB+LSB_f.*b_odd_sb;
                             
                             
                             
%                          
                      end 
                        
                        
                    %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                   
                    
      else
           ODD_LSB=zeros(1,odd_stop-odd_start+1);
          
          
      end     %����������ã�LSB��IF���
                        
                     
                    if sbnum_LSB<=sbants/2
                    sbnum_USB=sbnum_LSB+sbants/2;
                    else
                       sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
                    end
               
                
      sb_valid=sb_basic_data(sbnum_USB,6);
      
      if sb_valid   %����������ã��������Ӧ���� ��USB���侵��**********************************************     
                    
                  sb_x=sb_basic_data(sbnum_USB,1);
                  sb_y=sb_basic_data(sbnum_USB,2);
                  sb_z=sb_basic_data(sbnum_USB,3);
                  sb_a=sb_basic_data(sbnum_USB,4);
                  sb_p=sb_basic_data(sbnum_USB,5);
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_USB;   %���̲������λ����k
%                     sb_elev_max_usb=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min_usb=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                   USB_f=sb_Loss*Usb_Amp*Lsb_Amp*CSB_A*sb_a*cos(USB_func(odd_start:odd_stop)+phase_error); %*cos(phase_error)-sin(USB_func)*sin(phase_error));
                   
                    %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                        
                        
                     if checkSB_OBS   %���ڵ�
                             
                                ODD_USB=zeros(1,odd_stop-odd_start+1);
                           
%                         
                          
                       else
                           ODD_USB=USB_f.*b_odd_sb;
                             
%                          
                      end 
                       
                       
                    %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
             
                
               
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ߴ�
                     
%                     if sbnum_LSB<=sbants/2
%                     sbnum_USB=sbnum_LSB+sbants/2;
%                     else
%                        sbnum_USB=sbnum_LSB-sbants/2;    %ȡ����ӦUSB���ߺ�
%                     end
                        A_x=sb_x;
                  A_y=sb_y;
                  A_z=sb_z;
                  ZZ1=sb_z;

                       tt=(0-(-sb_z))/(fly_z-(-sb_z));
                    xx=sb_x+tt*(fly_x-sb_x);
                    yy=sb_y+tt*(fly_y-sb_y);
                    rr=sqrt(xx^2+yy^2);
                      
                      
               if rr<=counterpoint_R   %sb_elev_max_lsb
                
                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
                         sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
                          sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
                          sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
                          sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
                          reflect_rate=reflect_rate_couterpoint;
               else
%                    if  fly_angle<sb_elev_min_usb
                        sb_x=sb_basic_data_image_2_ground(sbnum_USB,1);
                          sb_y=sb_basic_data_image_2_ground(sbnum_USB,2);
                          sb_z=sb_basic_data_image_2_ground(sbnum_USB,3);
                          sb_a=sb_basic_data_image_2_ground(sbnum_USB,4);
                          sb_p=sb_basic_data_image_2_ground(sbnum_USB,5);
                          
                          
                        ZZ2=sb_z;
                      
                       %������㷨�Ǽ��ɻ��Ƿ����յ���ط�����źţ�������Ҫ�������㣬�������ߵ��ź��ܷ��ƹ��������䵽����ϡ�
                  %��һ���� ���㾵��������ɻ������ڴ���ϵ�����P
                  
                        P_zz=-((ZZ1-ZZ2)/2-ZZ1);
                        
                       tt=(P_zz-sb_z)/(fly_z-sb_z);
                   P_xx=sb_x+tt*(fly_x-sb_x);
                  P_yy=sb_y+tt*(fly_y-sb_y);
                    
               
                   
                 % �ڶ���������P����A��������ڷ������ϵĽ������꼰�뾶�����С�ڵ����뾶����ɻ����ղ������淴���źź͵��������ź�
%                   A_x=sb_x;
%                    A_y=sb_y;
%                  A_z=sb_z;
%                     ZZ1=sb_z;

                     tt=(0-A_z)/(P_zz-A_z);
                    xx=A_x+tt*(P_xx-A_x);
                    yy=A_y+tt*(P_yy-A_y);
                    rr=sqrt(xx^2+yy^2);
                       

                   if rr<=counterpoint_R   %sb_elev_max_lsb
                        
                       reflect_rate=0;
                   else
                      reflect_rate=reflect_rate_ground;
                   end  
                          
                      
%                    else   %���������MAX-MIN֮�䣬�����Ϊ���ӣ����о���������ʱ�ô���MAX��ֵ���㡣
%                        sb_x=sb_basic_data_image_1_counterpoint(sbnum_USB,1);
%                          sb_y=sb_basic_data_image_1_counterpoint(sbnum_USB,2);
%                           sb_z=sb_basic_data_image_1_counterpoint(sbnum_USB,3);
%                           sb_a=sb_basic_data_image_1_counterpoint(sbnum_USB,4);
%                           sb_p=sb_basic_data_image_1_counterpoint(sbnum_USB,5);
%                        
%                    end
               end 
                 d_sbant=sqrt((fly_x-sb_x)^2+(fly_y-sb_y)^2+(fly_z-sb_z)^2); 
                    loss_sb=32.45+20*log10(f_Usb1)+20*log10(d_sbant)-0-2.16-120-60;
                 sb_Loss=power(10,-loss_sb/20);
                 
                    delta_phase=(d_sbant-d_csb)*k_USB;   %���̲������λ����k
%                     sb_elev_max=atan(sb_z/(counterpoint_R-antenna_D/2))*180/pi;
%                     sb_elev_min=atan(sb_z/(counterpoint_R+antenna_D/2))*180/pi;
%                     b_even= rectpuls(t-w/2-(fix(i/2)-1)*w,w);
                    %    ֱ���ź� 
                    phase_error=Lsb_phase*pi/180+sb_p*pi/180+delta_phase;   %��λΪ����
                   USB_f=sb_Loss*Usb_Amp*Lsb_Amp*reflect_rate*CSB_A*sb_a*cos(USB_func(odd_start:odd_stop)+phase_error); %*cos(phase_error)-sin(USB_func)*sin(phase_error));
              %%%%%%%%      �ڵ�����  ������ڵ�������Ӧ�ź�Ϊ��       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %  function InObstacle=Inplane(p1,p2,Pttable)    %����boolֵ  ��2210��
                        p1=[fly_x,fly_y,fly_z];
                        p2=[sb_x,sb_y,sb_z];
                        checkSB_OBS=Inplane(p1,p2,Ptable,handles);
                     
                        if checkSB_OBS   %���ڵ�
                            
                            ODD_USB=ODD_USB+zeros(1,odd_stop-odd_start+1);
%                         
                          
                         else
                              ODD_USB=ODD_USB+USB_f.*b_odd_sb;
%                          
                       end 
                        
                    %%%%%%%%      �������         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %�д�����
                     %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
                   
                  
      else
             ODD_USB=zeros(1,odd_stop-odd_start+1);
         
      end    %�����������
       
%          LSB_ANT(i,:)=LSB.*b_odd.*b_cos;
%          USB_ANT(i,:)=USB.*b_odd.*b_cos;
    if i==1
        ODD_ANT=ODD_LSB+ODD_USB;
    else
        ODD_ANT=[ODD_ANT ODD_LSB+ODD_USB];
    end

 
  end   
    
     PerStr=fix(i/steps);
    waitstr=['processing.....',num2str(PerStr),'%','�� ',num2str(sbnum_LSB),' �űߴ�����'];
    waitbar(i/sbs_total,hwait,waitstr);
    
   
 end   
 


% --- Executes on button press in pushbutton_Theory.
function pushbutton_Theory_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Theory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
format long;
try
    % Multi-line editboxes are contained within a scroll-panel
      freq = get(handles.txt_freq,'String');  %ȡ��Ƶ��ֵ
    fmindex=get(handles.edit_FMI,'String');  %ȡ�õ�Ƶָ��
          d=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
        h_d=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
        D_r=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_c=get(handles.txt_Counterpoint_ReflectFactor,'String');  %����������txt_Counterpoint_ReflectFactor
  reflect_g=get(handles.txt_Ground_ReflectFactor,'String');  %����������
      csb_h=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       sb_h=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
      fly_r=get(handles.edit_R,'String');    %����뾶
        s_h=get(handles.edit_S_H,'String');    %��ʼ�߶�
        e_h=get(handles.edit_E_H,'String');    %�����߶�
        s_r=get(handles.edit_S_R,'String');    %��ʼ���롢��ʼ��λ
        e_r=get(handles.edit_E_R,'String');    %�������롢������λ
   fly_dial=get(handles.edit_Radial,'String');    %���澶��Ƕ�
   fly_step=get(handles.edit_step,'String');
   ph_ref=get(handles.edit_AZ_Align,'String');  %30Hz���࣬���ڷ�λУ��
   fm30_ph=get(handles.edit_FM30_Ph,'String');
   
    antenna_D=str2double( d);  %=get(handles.edit_AntennaArray_Dimension,'String');  %ȡ��������ֱ��
    couterpoint_H=str2double(h_d); %=get(handles.txt_Counterpoint_H,'String');     %ȡ�õ����߶�
    counterpoint_R=str2double(D_r); %=get(handles.txt_Counterpoint_R,'String');    %ȡ�õ����뾶
  reflect_rate_couterpoint=str2double(reflect_c); %=get(handles.txt_Couterpoint_Reflection,'String');  %����������
  reflect_rate_ground=str2double(reflect_g); %=get(handles.txt_Ground_Reflection,'String');  %����������
   CSB_H=str2double(csb_h); %=get(handles.CSB_H,'String');     %�ز����߸˸߶ȡ�
       SB_H=str2double(sb_h); %=get(handles.SBs_H,'String');     %�ߴ����߸˸߶ȡ�
     FlySimulate_Circle=str2double( fly_r); %=get(handles.edit_R,'String');    %����뾶
     start_H=str2double(   s_h);  %=get(handles.edit_S_H,'String');    %��ʼ�߶�
     end_H=str2double(   e_h);   %=get(handles.edit_E_H,'String');    %�����߶�
       start_range=str2double( s_r); %=get(handles.edit_S_R,'String');    %��ʼ���� ����λ
       end_range=str2double(  e_r); %=get(handles.edit_E_H,'String');    %�������롢��λ
   FlySimulate_Radial=str2double(fly_dial); %=get(handles.edit_Radial,'String');    %���澶��Ƕ�
    simulate_step=str2double(fly_step);  %���沽����Բ�ܵĲ�����λ��0.1�㣬������в�����λ��1�ס�
    az_align=str2double(ph_ref);
    Fm30_Ph=str2double(fm30_ph);
    freq_value=str2double(freq);
    fmi=str2double(fmindex);
 %%   
    if ~isnan(freq_value)
    wave_L=300/freq_value;
    half_wave_L=wave_L/2;
    array_d=fmi*wave_L/pi;
    else
          set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String','Ƶ����ֵ����');
    end
    
   sb=get(handles.radiobutton_SB48,'Value');
    if(sb==1)
        sbants=48;
    else
        sbants=50;
    end
    
 %%   

 
 
 
  
  T_on=1/30/sbants;  %ÿ���ߴ����߿���ʱ��Ƭ��     
 
  fc1= freq_value*1e6;  %�ز�Ƶ�ʣ���λHz
  fc=freq_value*1e6;
  f_Lsb=fc-9960;       %LSBƵ��
  f_Usb=fc+9960;       %USBƵ��
  
   f_Lsb1=fc1-9960;       %LSBƵ��
  f_Usb1=fc1+9960; 
  
  omega=30;            % �ز�����Ƶ��
  sucarrier=9960;      %���ز�Ƶ�� 9960
  
  ph30=az_align;              % ��Ƶ��λ
  ph30_FM=Fm30_Ph;
  
  
  ph_fc=0;             %�ز���λ
  R=50;                %�����迹
  
  T=0.1;   % ���γ���ʱ�䣬��λ����
                               %�������� fs=2^28;       %���ɵ����źŵ�ʱ������ʣ�
 fs=4e6;   %����Ƶ��,�������ԣ���4e6��6��8��10���Եõ��ԳƵ��ز�AM�ź� ����3��5e6��7��9���Գ�
% fs=500e3;
 N=fs*T;    %��������   N/T; %2^nextpow2(1*f);
  
 freq_rev=1/T;   %Ƶ�ʷֱ���=1/T=fs/N;
 

%  noisy=randn(1,N+1);
  dt=1/fs;
  t=0:dt:T-1/fs;    %����2��ʱ��Ƭ���ڵ��ź�;
  csb_power=str2double(get(handles.edit_CSB_Power,'String'));     % CSB����,��λW
  Lsb_phase=str2double(get(handles.edit_LSB_Phase,'String'));     %LSB��λ����λ�ȡ�
  Lsb_Amp=str2double(get(handles.edit_LSB_AMP,'String'))/100;     %LSB����
  Usb_Amp=str2double(get(handles.edit_LSB_USB_Ratio,'String'))/100;   % USB����
  AM30=str2double(get(handles.edit_AM30,'String'))/100;      % 30Hz AM ���ƶ�
  CSB_A=sqrt(2*csb_power*R);      
  
  csb_sideband=AM30*cos(2*pi*omega*t+ph30*pi/180);
% clr_sideband_150=n*sin(2*pi*omega2*t);
   
%    FFT_30=sin(2*pi*omega*t);
   
csb_mod=csb_sideband;
  
  CSB= CSB_A*(1+csb_mod).*cos(2*pi*fc*t+ph_fc*pi/180);
  LSB=Lsb_Amp*CSB_A*cos(2*pi*f_Lsb*t+Lsb_phase*pi/180);
  USB=Usb_Amp*Lsb_Amp*CSB_A*cos(2*pi*f_Usb*t+Lsb_phase*pi/180);
  SB_signal=LSB+USB;
%   sum_s=CSB+SB_signal;    %���ź�=�ز�+�ߴ�
  
  %����FM 9960 �ź�
   
  mtFM30= 1*cos(2*pi*omega*t+ph30_FM*pi/180);
  j_mt(1)=0;
  for i=1:length(t)-1   %����,Ƶ�ʵĻ�������λ��
      j_mt(i+1)=j_mt(i)+mtFM30(i)*dt;
  end
%   SUBCarrier_FM=AM30*cos(2*pi*sucarrier*t+2*pi*fmi*omega*j_mt);
  SUBCarrier_FM=Lsb_Amp*(1+Usb_Amp)*cos(2*pi*sucarrier*t+fmi*sin(2*pi*omega*t+ph30_FM*pi/180));
  Full_DVOR=CSB_A*(1+csb_mod+SUBCarrier_FM).*cos(2*pi*fc*t+ph_fc*pi/180);
 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����


xxlable="ʱ��";        %get(handles.edit_XLable,'String');
yylable="����";         %get(handles.edit_YLable,'String');
  figure(3);
  
  subplot(3,2,1);
plot(t,CSB);
% xlabel(xxlable);2
ylabel(yylable);
title(ftitle);
  t_p=T;
 axis([0 t_p -CSB_A*(1+AM30) CSB_A*(1+AM30)]);
 
   subplot(3,2,2);
    ftitle="�ߴ�";
plot(t,SB_signal);
 xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -CSB_A CSB_A]);
   

   
  
  
  
  v=get(handles.checkbox_Noise,'Value');
  
 for kk=1:200 
  
  if v==1
%   AM_signal=sum+randn(length(sum),1);%CSB+LSB+USB+noisy;
%   AM_signal=awgn(sum_s,20,'measured');   %��Ӹ�˹�������������źŵ�20dB
   AM_signal=awgn(Full_DVOR,18,'measured'); 
  else
%    AM_signal=sum_s;
   AM_signal=Full_DVOR;
   
  end
  
  
%    subplot(3,2,3);
%     ftitle="�ռ�������ź�";
% plot(t,AM_signal);
% % xlabel(xxlable);
% ylabel(yylable);
% title(ftitle);
%     t_p=T;
%  axis([0 t_p -2*CSB_A 2*CSB_A]);
  
  
  

  h=hilbert(AM_signal);   %����Hilbert�任������90��
%  yi=imag(h);       %Hilbert�任֮��õ��ĸ������鲿����Hilbert�任
%   xi=real(h);
 
  
  
  
 am_env1=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 
%  [y]=bandp(signal,hbandwidth,freq,srate,ntaps)
      % hbandwidth :��ͨ�˲����İ���������1KHz,������Ϊ500Hz
      % freq������Ƶ�� 
      %srate ����������
      %ntaps: filter coefficients. ���߽���taps. ��matlab�У����ǽ��� b
%    hbandwidth=5; %half bandwidth
%  function [y]=lowp(signal,lpf,srate,ntaps)
       % lpf:ͨ��Ƶ�� 
      %srate ����������
      %ntaps: filter coefficients. ���߽���taps. ��matlab�У����ǽ��� b
      % Butterworth IIR Filter
%   am_env=am_env1;
% %%
% am_env1=resample(am_env1,400e3,400e6);
%   am_env=filter(mylowpass,am_env1);
%   am_env=smooth(am_env1,100);
  am_env=am_env1;
%   
%  %clear AM_signal  SB_signal CSB USB LSB csb_mod j_mt Full_DVOR h SUBCarrier_FM  csb_sideband mtFM30;
%   
%   
%   fs=400e3;
%   dt=1/fs;
%   t=0:dt:T-1/fs;    %����2��ʱ��Ƭ���ڵ��ź�;
%   am_env=am_env(7000:35000);
%   DC=mean(am_env)
%   env_30AM=am_env-DC;
%   AM30_MOD=(max(env_30AM)-min(env_30AM))/DC/2; %FH(index_30)/FH(1);
%   
%   
%   
%   sucarrier=filter(mybandpass,am_env1);
%   env_9960= sucarrier(7000:end);
%     AM9960_MOD=(max( env_9960)-min( env_9960))/DC/2;
%     
%     diff9960=zeros(1,length(sucarrier));
%                                      for i=1:length(sucarrier)-1
%                                          diff9960(i)=(sucarrier(i+1)-sucarrier(i))/(1/fs);
%                                      end
%                                      
%                                      sss=abs(hilbert(diff9960)); 
%                                      sss=sss-mean(sss);
%      fm30=filter(mylowpass,sss);
%      
%     fm30=fm30(7000:35000);
% %     for i=1:1:length(fm30)
% %         if fm30(i)<0
% %             fm30(i)=0;
% %         end
% %         if fm30(i)>0
% %             fm30(i)=1;
% %         end
% %            if env_30AM(i)<0
% %             env_30AM(i)=0;
% %            end
% %         if env_30AM(i)>0
% %             env_30AM(i)=1;
% %         end
% %     end
%    
%                                      c_start=1;
%                                    c_stop=length(fm30)-c_start+1;  
%                                      maxId=1;
%                                    ww=blackman(c_stop-c_start+1);
% %                                  ww=blackmanharris(c_stop-c_start+1);
%                                      am30FFT=fft(env_30AM(c_start:c_stop).*ww');
% %                                     am30FFT=fft(am30_env(c_start:c_stop));
%                                    am30FFTAMP=abs(am30FFT);
%                                     [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
%                                     id30=id30+maxId+(30-10)/freq_rev-1;
%                                     ph30AM= angle(am30FFT(id30))*180/pi;
%                                     
%                                      fm30FFT=fft(fm30(c_start:c_stop).*ww');    %%%���ò�������������
% %                                         fm30FFT=fft(fm30_env(c_start:c_stop));
%                                    fm30FFTAMP=abs(fm30FFT);
%                                     [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
%                                    id30=id30+maxId+(30-10)/freq_rev-1;
%                                     ph30FM= angle(fm30FFT(id30))*180/pi;
%                                    
% 
% %                                      R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
% %                                      [Rmax,Rloc]=max(R);
% %                                     Rloc=Rloc-(c_stop-c_start+1);
% %                                     deg=Rloc*360*30/1250;
% 
%                                       deg=ph30FM-ph30AM
%    pp(kk)=deg;
%  end
%  figure(10);
%  plot(pp);
 %%
    %   ��������׼��λ��  -3.802538670044072
    %   �˲�������ʱ����ͬ�˲����Բ�ͬ�źŵ���ʱ��һ���������Ҫ����У����
%     
%     R=xcorr(fm30,env_30AM);
% [Rmax,Rloc]=max(R);
% Rloc=Rloc-length(fm30);
%   deg=Rloc*360*30/fs;
%   figure(10);
%   plot(fm30);
%   hold on;
%   plot(env_30AM);
    
%  am_env=lowp(am_env1,12000,fs,3);
% % % % % % % % % % % % % % % % % % %  [y]=butter_LP(x,
% wp1,ws1,r_p,r_s,fs)

%  am_env=butter_LP(am_env1,500,900,1,50,fs);
 
%  am_env=AM_signal;
           %  reshape(z,1,10000);
%  am_env=zz;

 %%%%%%%%%%%%-----------------------����FFT���ܣ�--------------------------
%    amp1=mean(am_env);
%   dft1=am_env.*FFT_30;   %�ֶ����㣬30Hz��DFT�任��ԭ������30Hz�źŶԲ����źŽ��о��������Ӧ����������
%   amp2=sum(dft1);        %���ݹ�ʽ���㣬
%   amp3=amp2/(N/2);       % �õ�30HzƵ�׳ɷݵķ��ȣ���FFT��ȣ�һ�¡���֤��ϡ�
  
%  subplot(4,1,3);
%   plot(t,zzz);
%   ylabel('����');xlabel('ʱ��');title('����첨');
 

 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 %fs=1/(T/(N-1));  %����Ƶ��
%  ffss=fs*linspace(0,1,NFFT);
 
 ft=fft(am_env,NFFT);
 
 FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 %fh1=FH.*conj(FH);;%conj()������y�����Ĺ������ʵ���Ĺ������������
% ff=fs*(0:NFFT/2-1);% F F T �任���Ӧ��Ƶ�ʵ�����
 %b=FH(90);
 
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
  
  index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 
 
  AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);

  
     subplot(3,2,3);
    ftitle="�ռ�������ź�";
plot(t,AM_signal);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
    t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);

   
    
   subplot(3,2,4);
    ftitle="AM�����  AM30="+num2str(AM30_MOD,'%1.4f')+"  AM9960="+num2str(AM9960_MOD,'%1.4f');
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
   t_p=T;
 axis([0 t_p -2*CSB_A 2*CSB_A]);
  
  str1={'ֱ��:', '30HzAM:' , '9960HzAM:'};
 t_p=0.07;
   
  
 
 text(t_p/2,-30, str1,'Color','red','FontSize',8);
 str2={num2str(FH(1)),  num2str(FH(index_30),'%1.4f'),  num2str(FH(index_9960),'%1.4f')};
  text(t_p/2+40*t_p/80,-30,str2,'Color','red','FontSize',8);
 
   subplot(3,2,5);
    ftitle="Ƶ��ͼ";
  freqaxis=(-NFFT/2:NFFT/2-1)*freq_rev;   %fshift = (-n/2:n/2-1)*(fs/n)
  YY=fftshift(FH);
  plot(freqaxis,YY);

xlabel("Ƶ��");
ylabel("����");
title(ftitle);




% 
%  h=hilbert(AM_signal);   %����Hilbert�任������90��
% 
%  am_env=abs(h);       %sqrt(yi.*yi+xi.*xi); Hilbert�任ʵ���ǰ���첨
 


 
 NFFT=N; %                    2^nextpow2(length(zzz));  %));%�ҳ�����y�ĸ���������2��ָ��ֵ FFT����
 ft=fft(am_env,NFFT);
  FH=abs(ft)*2/NFFT; %��f�źŽ���DFT���õ�Ƶ�ʵķ�ֵ�ֲ�
 FH(1)=FH(1)/2;   % DCֱ���ɷ�
 index_30=30/freq_rev+1;
 index_9960=9960/freq_rev+1;
 AM30_MOD=FH(index_30)/FH(1);
 AM9960_MOD=FH(index_9960)/FH(1);
 ang30=angle(ft(index_30))*180/pi;
 
 
 
   y_data_dbm = 10*log10((FH.^2)/50/2)+30;  %���㹦��dBmֵ,���ݷ��ȼ��㹦�ʣ���1/2 A^2/R,��λ��dBW,����30������dBm��
  [~,maxId]=max(FH(1:NFFT));
    rflevel=y_data_dbm(maxId);
     RF_Level= num2str(rflevel);  
      
     [~,id30]=max(FH(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
     
  [~,id9960]=max(FH(maxId+(9960-900)/freq_rev:maxId+(9960+900))/freq_rev);
   a9960=FH(maxId+(9960-900)/freq_rev-1+id9960);
   AM9960=FH(maxId+(9960-900)/freq_rev+id9960-1)/FH(1)/(50/30);                               
                                     
                                     maxSuId_1=0;
                                      maxSuId_2=0;
                                      
                                      
                                    for iid=maxId+(9960-1000)/freq_rev:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_1=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    for iid=maxId+(9960+1000)/freq_rev:-1:maxId+9960/freq_rev
                                        if FH(iid)>(a9960/2)
                                         maxSuId_2=iid;
                                         break;
                                        end
                                        
                                    end
                                    
                                    ss=ifft(ft(maxId),NFFT);
                                   sss=real(ss);
                                   DC1=mean(sss);
                                   
                                   
                                    
                                     ss=ifft(ft(maxId:maxId+(30+210)/freq_rev),NFFT);
                                     am30_env=real(ss);
                                     am30_env=am30_env-mean(am30_env);
                                      
                                     AM30=2*max(am30_env(50000/freq_rev:NFFT-50000/freq_rev))/DC1;%FindAmp(30,nfft,sss);
                                      
                                   
% F   =  [0:0.05:0.95]; 
% A  =  [1    1      0     0     0    0      0     0     0    0     0     0     0     0     0     0    0   0   0   0] ;
% b  =  firls(20,F,A);
% 
%                                      Signal_Filter= filter(b,a,am30_env);
%                                      am30_env=Signal_Filter;
%                                      
%                                   
%                                    
%                                    am30_env=am30_env.*1000/5;
                                   
%                                      disp(['30HzAM: ',num2str(A_30/DC1)]);
                                     
%                                     sss=resample(yy,2^15,2^18);
%                                      fm30=lowp(abs(sss),40,90,1,30,2^15);%fft(sss,nfft);
%                                     
%                                      
%                                    am30_env=abs(fm30);
                                     ft(maxId+(30-20)/freq_rev:maxId+(30+100)/freq_rev)=0; %(ft(maxId+(30+60)/freq_rev-1)+ft(maxId+(30+60)/freq_rev+1))/2;

                                      ss=ifft(ft(maxId+(9960-600)/freq_rev:maxId+(9960+600)/freq_rev),NFFT);
%                                       ss=ifft(ft,NFFT);
                                   am9960_env=real(ss);
                             
                                   AM9960=2*max(am9960_env(5000/freq_rev:NFFT-5000/freq_rev))/DC1;%FindAmp(9960,nfft,sss);
%                                      disp(['9960HzAM: ',num2str(A_9960/DC1)]);
                                     
                                     
                                        diff9960=zeros(1,length(am9960_env));
                                     for i=1:length(am9960_env)-1
                                         diff9960(i)=(am9960_env(i+1)-am9960_env(i))/(1/fs);
                                     end
                                     
                                     sss=abs(hilbert(diff9960));


% sss=fmdemod(am9960_env,9960,rtlsdr_fs,480);    %diff(am9960_env);
                                    
%                                      sss=resample(rtl_fft,2^15,2^18);
%                                      fm30=lowp(sss,40,90,1,30,2^15);%fft(sss,nfft);
                                     fm30_env=real(sss)-mean(real(sss));
%                                      figure(10);
%                                      plot(filtfilt(b,a,fm30_env));
%                                      fm30_env(1)=fm30_env(2);
%                                      am30_env=resample(am30_env,2^15,nfft);
%                                      fm30_env=resample(fm30_env,2^15,nfft);
%                                      Wc=2*50/new_nfft;                                          %��ֹƵ�� 50Hz
%                                      [b,a]=butter(4,Wc);
%                                      Signal_Filter=filter(b,a,fm30_env);
%                                      
%                                      fm30_env=Signal_Filter;
                                     
                                        am30_env=resample(am30_env,1250,NFFT);
                                     fm30_env=resample(fm30_env,1250,NFFT);
                                     
                                       c_start=1;
                                   c_stop=length(am30_env)-c_start+1;  
                                     
                                   ww=blackman(c_stop-c_start+1);
%                                  ww=blackmanharris(c_stop-c_start+1);
                                     am30FFT=fft(am30_env(c_start:c_stop).*ww');
%                                     am30FFT=fft(am30_env(c_start:c_stop));
                                   am30FFTAMP=abs(am30FFT);
                                    [~,id30]=max(am30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                    id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30AM= angle(am30FFT(id30))*180/pi;
                                    
                                     fm30FFT=fft(fm30_env(c_start:c_stop).*ww');    %%%���ò�������������
%                                         fm30FFT=fft(fm30_env(c_start:c_stop));
                                   fm30FFTAMP=abs(fm30FFT);
                                    [~,id30]=max(fm30FFTAMP(maxId+(30-10)/freq_rev:maxId+(30+10)/freq_rev));
                                   id30=id30+maxId+(30-10)/freq_rev-1;
                                    ph30FM= angle(fm30FFT(id30))*180/pi;
                                   

%                                      R=xcorr(am30_env(c_start:c_stop),fm30_env(c_start:c_stop));
%                                      [Rmax,Rloc]=max(R);
%                                     Rloc=Rloc-(c_stop-c_start+1);
%                                     deg=Rloc*360*30/1250;

                                      deg=ph30FM-ph30AM;
                                      
                                   
            
                                    if deg<0
                                        deg=deg+360;
                                    end
%                                     if deg<-180
%                                         deg=deg+360;
%                                     end
%                                     
%                                     ang_r=ang;
%                                     if ang>180
%                                         ang_r=ang-360;
%                                     end
                                    ang=ph30_FM-ph30; %Fm30_Ph;
                                    az_error=deg-ang;  %���㷽λ���
                                    
%                                     az_error=az_error-180; 
                                    if az_error>180
                                        az_error=az_error-360;
                                    end
                                    if az_error<-180
                                        az_error=az_error+360;
                                    end
                                    
%                                         if abs(az_error)>=180
%                                             az_error=az_error-360;
%                                         end
                                    
%                                     figure(2);
%                                     subplot(311);
%                                     plot(am30_env);
%                                    title('30Hz AM�ź�');
%                                    subplot(312);
%                                    plot(fm30_env);
%                                     title('30Hz FM�ź�');
%                                     subplot(313);
%                                     plot(R);
%                                     title(['XCORR ','Max Rloc==',num2str(Rloc)]);
 
                                      fmi=(maxSuId_2- maxSuId_1)*freq_rev/2.0/30.00;
                                     
                                      

                                          
                                          VOR_AZ=num2str(az_error,'%.4f');
                                          VOR_30HzAM=num2str(round(AM30*100*1000)/1000);
                                          VOR_9960HzAM=num2str(round(AM9960*100*1000)/1000);
                                          VOR_FMI=num2str(fmi);
                                     results=[ang,rflevel,az_error,round(AM30*100*1000)/1000,round(AM9960*100*1000)/1000,fmi];
     
     figure(3);
   subplot(3,2,6);
%     ftitle="AM������-AM30MOD="+VOR_30HzAM+"  AM9960MOD="+ VOR_9960HzAM;
      ftitle="AZ: "+num2str(deg,'%.3f')+" | AM30="+VOR_30HzAM+" | AM9960="+ VOR_9960HzAM+" | AZ_E_r_r="+VOR_AZ+'  FMI:'+VOR_FMI;
plot(t,am_env);
xlabel(xxlable);
ylabel(yylable);
title(ftitle);
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     pp(kk)=az_error;
     
     
     
 end
 figure(10);
 plot(pp);
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ftitle="�ز�"; %get(handles.edit_figureTitle, 'String');  %ͼ�α�����
% xxlable="ʱ��";        %get(handles.edit_XLable,'String');
% yylable="����";         %get(handles.edit_YLable,'String');
%   figure(3);
  
%  
%    
  


 







catch ErrorInfo
    % probably a single-line editbox
  %  throw(ErrorInfo);  %��ֹ����ִ��
%   disp(ErrorInfo);
    
     set(handles.txt_Error,'Visible','On');
    set(handles.txt_Error,'String',ErrorInfo.message+"�к�:"+string(ErrorInfo.stack(1).line)+"��������"+ErrorInfo.stack(1).name);
end



function edit_FM30_Ph_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FM30_Ph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FM30_Ph as text
%        str2double(get(hObject,'String')) returns contents of edit_FM30_Ph as a double


% --- Executes during object creation, after setting all properties.
function edit_FM30_Ph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FM30_Ph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%����Ӧ�˲���LMS
%��С�����㷨��Least Mean Square��LMS���� Bernard Widrow �� Marcian E. Hoff ��������������˲�������ʹ�����Mean Square Error��MSE���ﵽ��С��
%LMS�㷨����Ϊ�ǻ���ѧϰ���������Ҳ�Ƚ����õ��㷨���������жԲ�����ѧϰʹ�õľ���LMS��˼�룬��ͨ���źŴ�������LMSҲ�ǳ���������������Ӧ�˲�����
 % �������:
%     xn   ������ź�����      (������)
%     dn   ����������Ӧ����    (������)
%     M    �˲����Ľ���        (����)
%     mu   ��������(����)      (����)     Ҫ�����0,С��xn����ؾ����������ֵ�ĵ���    
% �������:
%     W    �˲�����Ȩֵ����     (����)
%          ��СΪM x itr,
%     en   �������(itr x 1)    (������)  
%     yn   ʵ���������         (������)
function [yn,W,en]=LMS(xn,dn,M,mu)
   itr = length(xn);
en = zeros(itr,1);             % �������,en(k)��ʾ��k�ε���ʱԤ�������ʵ����������
W  = zeros(M,itr);             % ÿһ�д���һ����Ȩ����,ÿһ�д���-�ε���,��ʼΪ0
% ��������
for k = M:itr                  % ��k�ε���
    x = xn(k:-1:k-M+1);        % �˲���M����ͷ������
    y = W(:,k-1).' * x;        % �˲��������
    en(k) = dn(k) - y ;        % ��k�ε��������
    % �˲���Ȩֵ����ĵ���ʽ
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end
% ������ʱ�˲������������  r���û��yn���ز������Բ�Ҫ�����
yn = inf * ones(size(xn)); % inf ����������˼
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;%�����õ�����ѹ��Ƶõ����
end

% �ź��˲�Ӧ�þ���
% % % % % % % % % % % % % % % % xn = xs+xn;
% % % % % % % % % % % % % % % % xn = xn.' ;   % �����ź�����
% % % % % % % % % % % % % % % % dn = xs.' ;   % Ԥ�ڽ������
% % % % % % % % % % % % % % % % M  = 20 ;   % �˲����Ľ���
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % rho_max = max(eig(xn*xn.'));   % �����ź���ؾ�����������ֵ
% % % % % % % % % % % % % % % % mu = (1/rho_max) ;    % �������� 0 < mu < 1/rho
% % % % % % % % % % % % % % % % [yn,W,en] = LMS(xn,dn,M,mu);


    function [y e]=my_LMS(x,d,u,N,e)   
      % [y e]=LMS(x,d,u,N,e) 
      %u=2*�������ӣ� һ����������Ϊ0.0022,��2e-3
      %w:���Ƶ�FIR�˲���
      %y:�������y[n]
      %x:��������x[n];
      %d:�������飬���ȵ���x[n]
      %u:���ڲ���
      %N:FIR�˲�������
      
        
        
        
   function [y]=bandp(signal,hbandwidth,freq,srate,ntaps)
      % hbandwidth :��ͨ�˲����İ���������1KHz,������Ϊ500Hz
      % freq������Ƶ�� 
      %srate ����������
      %ntaps: filter coefficients. ���߽���taps. ��matlab�У����ǽ��� b
%    hbandwidth=5; %half bandwidth
nyquist=srate/2;
% ntaps=200;
hpf=(freq-hbandwidth)/nyquist;
lpf=(freq+hbandwidth)/nyquist;
[b,a]=fir1(ntaps, [hpf lpf]); %calc filter coefficients
signal_f=filter(b,a,signal); %filter signal
 y=signal_f;
 
 
 function [y]=lowp(signal,lpf,srate,ntaps)
       % lpf:ͨ��Ƶ�� 
      %srate ����������
      %ntaps: filter coefficients. ���߽���taps. ��matlab�У����ǽ��� b
      % Butterworth IIR Filter
% srate=1000; % sample rate
% npts=2000; %npts in signal
Nyquist=srate/2; %Nyquist frequency
% lpf=300; %low-pass frequency
order=ntaps; %filter order
% t=[0:npts-1]/srate; %time scale for plot
% x=(rand(npts,1)*2)-1; % raw data from -1 to +1

[b,a]=butter(order,lpf/Nyquist); %create filter coefficients

y=filter(b,a,signal); % filter using 'b' and 'a' coefficients

     function   [y]=butter_LP(x, wp1,ws1,r_p,r_s,fs)
%���һ�������ַ��ͨ�˲���,Ҫ���50Hz��Ƶ�ʷ�������,���������˵�
wp = wp1/(fs/2);  %ͨ����ֹƵ��,ȡ50~100�м��ֵ,�������һ��(65Hz)
ws = ws1/(fs/2);  %�����ֹƵ��,ȡ50~100�м��ֵ,�������һ��(85Hz)
alpha_p = r_p; %ͨ���������˥��Ϊ  db
alpha_s = r_s;%���������С˥��Ϊ  db
%��ȡ�����ͽ�ֹƵ��
[ N1,wc1 ] = buttord( wp , ws , alpha_p , alpha_s);
%���ת�ƺ���ϵ��
[ b,a ] = butter(N1,wc1,'low');
%�˲�
y = filter(b,a,x);
 
     function   [y]=butter_HP(x, wp1,ws1,r_p,r_s,fs)
%���һ����ͨ�˲���,Ҫ���400Hz��Ƶ�ʷ�������,���������˵�
wp = wp1/(fs/2);  %ͨ����ֹƵ��,ȡ200~400�м��ֵ,�������һ��
ws = ws1/(fs/2);  %�����ֹƵ��,ȡ200~400�м��ֵ,�������һ��
alpha_p = r_p; %ͨ���������˥��Ϊ  db
alpha_s = r_s;%���������С˥��Ϊ  db
%��ȡ�����ͽ�ֹƵ��
[ N2,wc2 ] = buttord( wp , ws , alpha_p , alpha_s);
%���ת�ƺ���ϵ��
[ b,a ] = butter(N2,wc2,'high');
%�˲�
y = filter(b,a,x);

     function   [y]=butter_BP(x, wp1,wp2,ws1,ws2,r_p,r_s,fs)
%���һ����ͨ�˲���,Ҫ���50Hz��400Hz��Ƶ�ʷ����˵�,������������
wp = [wp1 wp2 ] / (fs/2);  %ͨ����ֹƵ��,50~100��200~400�м��ȡһ��ֵ,�������һ�� 65-385Hz
ws = [ws1 ws2 ] / (fs/2);  %�����ֹƵ��,50~100��200~400�м��ȡһ��ֵ,�������һ�� 75-375Hz
alpha_p = r_p; %ͨ���������˥��Ϊ  db
alpha_s = r_s;%���������С˥��Ϊ  db
%��ȡ�����ͽ�ֹƵ��
[ N3,wn ] = buttord( wp , ws , alpha_p , alpha_s);
%���ת�ƺ���ϵ��
[ b,a ] = butter(N3,wn,'bandpass');
%�˲�
y = filter(b,a,x);

     function   [y]=butter_BS(x, wp1,wp2,ws1,ws2,r_p,r_s,fs)
%���һ�������˲���,Ҫ���50Hz��400Hz��Ƶ�ʷ�������,���������˵�
wp = [wp1 wp2 ] / (fs/2);  %ͨ����ֹƵ��?,50~100��200~400�м��ȡһ��ֵ,�������һ��
ws = [ws1 ws2 ] / (fs/2);  %�����ֹƵ��?,50~100��200~400�м��ȡһ��ֵ,�������һ��
alpha_p = r_p; %ͨ���������˥��Ϊ  db
alpha_s = r_s;%���������С˥��Ϊ  db
%��ȡ�����ͽ�ֹƵ��
[ N4,wn ] = buttord( wp , ws , alpha_p , alpha_s);
%���ת�ƺ���ϵ��
[ b,a ] = butter(N4,wn,'stop');
%�˲�
y = filter(b,a,x);
         
         

% --- Executes during object creation, after setting all properties.
function pushbutton_Theory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_Theory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function [max_f, min_f, mean_f,max_indix,min_indix]=frequencycounter(in_data,in_fs)
%           pulse_peak=max(yy_data);
p=20e6;
q=in_fs;

yy_data=resample(in_data,p,q);
fs=in_fs*p/q;


          threshold=0;

%           rise1=0.1*pulse_peak;
%           rise2=0.9*pulse_peak;
% if yy_data(1)<threshold    
% flag=0;
%  
% else    
% flag=1;
%  
% end
flag=0;

old=yy_data(1);
% new=yy_data(2);
 
max_t=0;
min_t=0;
mean_t=0;
max_indix=0;
min_indix=0;
% t1=0;
% t2=0;

for t=2:1:length(yy_data)  %   
new=yy_data(t);   

 %������  
 
        if(new>threshold)&&(old<=threshold)&&(flag==0)        
            flag=1;  
            t1=t-abs(new)/(new-old);
        end
        
        if (new>threshold)&&(old<=threshold)&&(flag==1)
                 t2=t-abs(new)/(new-old);
                f1=t2-t1;
                t1=t2;
                if max_t<f1
                    max_t=f1;
                    max_indix=t;
                end
                if min_t==0 || min_t>f1
                        min_t=f1;
                        min_indix=t;
                end
                if mean_t==0
                    mean_t=f1;
                else
                    mean_t=(mean_t+f1)/2;
                end
        end
 




old=new;

end
max_f=fs/min_t;
min_f=fs/max_t;
mean_f=fs/mean_t;
   
