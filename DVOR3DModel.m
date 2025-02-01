function varargout = DVOR3DModel(varargin)
% DVOR3DModel MATLAB code for DVOR3DModel.fig
%      DVOR3DModel, by itself, creates a new DVOR3DModel or raises the existing
%      singleton*.
%
%      H = DVOR3DModel returns the handle to a new DVOR3DModel or the handle to
%      the existing singleton*.
%
%      DVOR3DModel('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DVOR3DModel.M with the given input arguments.
%
%      DVOR3DModel('Property','Value',...) creates a new DVOR3DModel or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DVOR3DModel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DVOR3DModel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DVOR3DModel

% Last Modified by GUIDE v2.5 15-Apr-2020 11:04:23

% Begin initialization code - DO NOT EDIT
clc;

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DVOR3DModel_OpeningFcn, ...
                   'gui_OutputFcn',  @DVOR3DModel_OutputFcn, ...
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

% --- Executes just before DVOR3DModel is made visible.
function DVOR3DModel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DVOR3DModel (see VARARGIN)

% Choose default command line output for DVOR3DModel
handles.output = hObject;
 
handles.ht=timer;  %���ʱ�ӣ���ʱ��������ˢ�»���
set(handles.ht,'ExecutionMode','FixedRate');
set(handles.ht,'Period',1);     % ÿ��ˢ��һ��
set(handles.ht,'TimerFcn',{@dispnow,handles});  %�󶨼�ʱ��ִ��FUNCTION
% Update handles structure
guidata(hObject, handles);
start(handles.ht);           %������ʱ��

handles.DVOR=varargin;  %����Ĳ�����  DVOR3DModel(CSB_H,antenna_D,counterpoint_R,couterpoint_H,sb_basic_data);  %�����ز����߷��ߣ�������ֱ���������뾶�������߶ȣ��ߴ����������

 guidata(hObject, handles)
 %%%%��ʼ��uiTable�е����ݣ���2��planeΪ��
 
  dt1=cell(2,7);
                        dt1{1,1}=500;   %����
                        dt1{1,2}=30;    %��λ
                        dt1{1,3}=40;    %�ϱ߿�
                        dt1{1,4}=60;    %�±߶�
                        dt1{1,5}=80;    %�߶ȣ���ԭʼ�߶ȵ��ڸ�ֵ����TAN(��ǣ�
                      
                        dt1{1,6}=5;     %%Z������ƫ�� %3.3f
                        dt1{1,7}=20;    %Z����ת
                        dt1{1,8}=3;     %��б�Ƕ�
                        dt1{1,9}=1.0;   %����ϵ��
                        dt1{1,10}=true; %���ñ��
       
                        dt1{2,1}=800;
                        dt1{2,2}=90;
                        dt1{2,3}=60;
                        dt1{2,4}=60;
                        dt1{2,5}=80;
                        dt1{2,6}=5; % %3.3f
                        dt1{2,7}=20;
                        dt1{2,8}=3;
                        dt1{2,9}=1.0;
                        dt1{2,10}=true;
 set(handles.uitable1,'Data',dt1);
 
draw_DVOR_BASIC(hObject,eventdata,handles);  %��DVOR����������λ��

% This sets up the initial plot - only do when we are invisible
% so window can get raised using DVOR3DModel.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes DVOR3DModel wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function dispnow(hObject,eventdata,handles)
set(handles.disptime,'String',datestr(now));
% draw_DVOR_BASIC(hObject,eventdata,handles);


%%  ��DVOR̨
function draw_DVOR_BASIC(hObject,eventdata,handles)

 %���ɱߴ���������  1�����ߺţ�2���Ƕȣ�3�����룬4���߶ȣ�5�����ȣ�6����λ��7������
   %�����ز����߷��ߣ�������ֱ���������뾶�������߶ȣ��ߴ����������

  %         1          2          3        4           5  (48*6)     
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% % %   DVOR3DModel(CSB_H,antenna_D,counterpoint_R,couterpoint_H,sb_basic_data,handles);  %�����ز����߷��ߣ�������ֱ���������뾶�������߶ȣ��ߴ��������,�����ھ����
  
if isempty( handles.DVOR)     
   
    return;
end
global figOption;
 global plotCount;
 global Ptable;
 
axes(handles.axes1);    %ȡ����ͼ������
cla;                   %�����ҵ��


% n=[1 1 1];   %;?%������n
% r=handles.DVOR{2};;  %?%Բ�İ뾶Ϊ1
% c=[1 1 1];   %?%Բ�ĵ�����
% theta=(0:2*pi/100:2*pi)'; %;?%theta�Ǵ�0��2*pi
% a=cross(n,[1 0 0]); %?%n��i��ˣ���ȡa����
% if ~any(a)  %?%���aΪ����������n��j���
%  a=cross(n,[0 1 0]);
% end
% b=cross(n,a); %��ȡb����
% a=a/norm(a); %��λ��a����
% b=b/norm(b); %��λ��b����
% 
% c1=c(1)*ones(size(theta,1),1);
% c2=c(2)*ones(size(theta,1),1);
% c3=c(3)*ones(size(theta,1),1);
% 
% x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%Բ�ϸ����x����
% y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%Բ�ϸ����y����
% z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%Բ�ϸ����z����
% 
% plot3(x,y,z)
% xlabel('x��')
% ylabel('y��')
% zlabel('z��')

h_parent=handles.DVOR{6};

N=200;                % ��Բ�ķ��ű��ʣ� 360��/200
r1=handles.DVOR{3};   %�������뾶
h=handles.DVOR{4};    %�������߶�
sbants=handles.DVOR{5};  %�ߴ�������������Ƕȣ��߶ȣ���������
sb=sbants(:,1:3);        %�ߴ����ߵ����� x,y,z

rr=linspace(0,r1,N);   %  ���򲽽�����
theta=linspace(0,2*pi,N); %�����Ƕ�
[R,the]=meshgrid(rr,theta);  %���ɾ����൱��N��ͬ��Բ�����ֳ�N������������R���󣬺�the����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %        R=[0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
% %           [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%           [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=0*R.^0;   %����Z��0����
[x,y]=pol2cart(the,R);   %������תֱ������
mesh(x,y,Z,'edgecolor',[0.2,0.5,1],'EdgeAlpha',0.5);    % �����磬��������Բ
 %,'FaceColor',[1 1 1]);
hold on;
% plot3(0,0,-2*rr ,'.');
% plot3(0,0,2*rr,'.')
colormap Colorcube;

%colorMapNames={'Parula','Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter','Gray','Bone','Copper','Pink','Lines','Colorube',
%                 'Prism','flag','white'

r2=handles.DVOR{2}/2;   %������ֱ�����뾶=D/2
rr=linspace(0,r2,N);    
theta=linspace(0,4*pi,N);  % 4pi,��2Ȧ
[R,the]=meshgrid(rr,theta);
Z=0.101*R.^0;          %���������Բ��Z������ƽ���ϣ���Ȼ��ɫ���ص�
[x,y]=pol2cart(the,R);   %������תֱ������
C=x.*y;
mesh(x,y,Z,'edgecolor','y','EdgeAlpha',0.1);  %,'FaceColor',[0,1,1]);
 colormap Spring;
 
 p1=[0,0,2*r1];
 p2=[0,0,-h];
 XX=[p1(1),p2(1)];
  YY=[p1(2),p2(2)];
   ZZ=[p1(3),p2(3)];
plot3(XX,YY,ZZ,...             %%%��һ��ͨ�������ز����ߣ���Բ�ĵĴ�ֱ�ߣ���Z��2*r1��-�裬ȡ�����档DVOR������ԭ�㶨�ڷ�����������ģ��������ز����ߵĻ����㡣
         'LineWidth',2,...
          'MarkerEdgeColor','r',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10);

[a,b]=size(sb);    %��48���ߴ����ߣ�����һ����0=����һ����ֱ�߱�ʾ��
for i=1:a
    plot3(sb(i,1),sb(i,2),sb(i,3),'o','Color',[0,1,0],'LineWidth',2);   %��һ����o'��״������ALFORD��������
    plot3([sb(i,1),sb(i,1)],[sb(i,2),sb(i,2)],[0,sb(i,3)],'LineWidth',1,'Color',[1 0 0]);  %����ֱ�ߣ�����ߴ����߸�
end


%%%%%%%%%%%%%���ɻ�λ��----------------------------------
pp=get(handles.checkbox_plane,'Value');
CSB_X=0;
CSB_Y=0;
CSB_Z=sb(1,3);

if pp
  fly_z=str2double(get(handles.edit_h,'String'))-h;   %�ɻ��߶ȼ�ȥ�����߶ȣ���Z��߶�
  FlySimulate_Circle=str2double(get(handles.edit_distance,'String'));
  ang=str2double(get(handles.edit_az,'String'));
  
    if fly_z<0
       
       set(handles.edit_cross,'String',"�ɻ��ڵ����·������������÷���ģʽ�Ͳ�����");   %�����ֵΪ�������˳���
    return;
    end
    
    fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
%     fly_angle=atan(fly_z/fly_d)*180/pi;   %����ɻ������ǣ���λ�Ƕȡ�
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
   plot3(fly_x,fly_y,fly_z,'o','Color',[1,0,0],'LineWidth',10);  %��һ���ű�ʾ�ɻ���
   plot3([0,fly_x],[0,fly_y],[CSB_Z,fly_z],'LineWidth',1,'Color',[1 0 0]);  %���ز�������ɻ������� 
end
%%%%%%%%%%%%------------����-���ɻ�λ��


% % % 
% % % %%%%%%%%%%%%%%%%��ȡ�ϰ�������������ϰ���ƽ�棬���λ򳤷��Σ������ϱ߿��±߿��߶ȣ���б����ת�Ȳ�������
% % %    dt1{1,1}=500;   %����
% % %                         dt1{1,2}=30;    %��λ
% % %                         dt1{1,3}=40;    %�ϱ߿�
% % %                         dt1{1,4}=60;    %�±߶�
% % %                         dt1{1,5}=80;    %�߶ȣ���ԭʼ�߶ȵ��ڸ�ֵ����TAN(��ǣ�
% % %                       
% % %                         dt1{1,6}=5;     %%Z������ƫ�� %3.3f
% % %                         dt1{1,7}=20;    %Z����ת
% % %                         dt1{1,8}=3;     %��б�Ƕ�
% % %                         dt1{1,9}=1.0;   %����ϵ��
% % %                         dt1{1,10}=true; %���ñ��

dt1=get(handles.uitable1,'Data');  %��ȡ�ϰ������е����ݱ�
[a,b]=size(dt1);
iii=1;
for i=1:a       %ѭ�������ϰ���ƽ�棬ÿ��һ�����õ�N���ϰ�������ƽ�ϵ����ꡣ����ֵ�ڹ�������Ptable.
                       dist=dt1{i,1};   %����,��ͶӰ��XOYƽ��ľ��룬б�����
                       az= dt1{i,2};    %��λ
                        top_width=dt1{i,3};    %�ϱ߿�
                        bottom_width=dt1{i,4};    %�±߶�
                        elev_h=dt1{i,5};    %�߶ȣ���ԭʼ�߶ȵ��ڸ�ֵ����TAN(��ǣ�
                      
                        Z_shift=dt1{i,6};     %%Z������ƫ��,��
                        Z_rotate=dt1{i,7};    %Z����ת
                        Z_inclination=dt1{i,8};     %Z�������б�Ƕ�
                        reflect_factor=dt1{i,9};   %����ϵ��
                        active=dt1{i,10};          %���ñ��
 
if active==true     %%�Ƚ�ƽ���Ƶ�������ԭ�㴦��������x=0��ƽ���ϣ�������Ϊy����,��Ϊx����,ֽ�ⷽ��ΪZ����
                    origin_H=elev_h/cos(Z_inclination*pi/180);     %  ����ƽ�治��бʱ�ĸ߶�
                    A0=[0;bottom_width/2;0];     %���������õ��Ƿֺţ��������е�����ʽ�� ��ʾX Y Z .���Σ����Σ����ĸ�������꣬A�����ϣ�B�����£�C����B�Ϸ���D����A�Ϸ�
                    B0=[0;-bottom_width/2;0];
                    C0=[0;-top_width/2;origin_H];
                    D0=[0;top_width/2;origin_H];
                    
                    
                    %  Z�᷽����б����ʱAB���䣬CD���⡢����б����A0��B0���䣬ʵ����Y����ת
                    shift_Y=[cos(Z_inclination*pi/180),0,sin(Z_inclination*pi/180);0,1,0;-sin(Z_inclination*pi/180),0,cos(Z_inclination*pi/180)];
                     
                    A0=shift_Y*A0;
                    B0=shift_Y*B0;
                    C0=shift_Y*C0;
                    D0=shift_Y*D0;
                    
                    
                    %��Z����ת
                    shift_Z=[cos(Z_rotate*pi/180),-sin(Z_rotate*pi/180),0;sin(Z_rotate*pi/180),cos(Z_rotate*pi/180),0;0,0,1];
                    
                    A0=shift_Z*A0;
                    B0=shift_Z*B0;
                    C0=shift_Z*C0;
                    D0=shift_Z*D0;
                    
                    
                    
                    %ƽ�Ƶ� dist,Az��
                    shift_T=[dist*sin(az*pi/180);dist*cos(az*pi/180);Z_shift];
                    shift_y1=[0;bottom_width/2;0];
                    shift_y2=[0;top_width/2;0];
                    
                    
                    
                    
                    
                    
                    A=A0+shift_T;
                    B=B0+shift_T;
                    C=C0+shift_T;
                    D=D0+shift_T;
                    
                    PointTable(iii,1)={[A,B,C,D]};   %��ӵ���������.����Ԫ������cell�������ĸ�����ͷ���ϵ����ֻ��ACTIVE���ϰ������ӣ�����ACTIVE�Ĳ����
                     PointTable(iii,2)={reflect_factor};
                     iii=iii+1;
                     
                    %                     h_parent.Ptable=PointTable;
%                     guidata(h_parent,h_parent.Ptable);
                    
%                     A=[100;0;0];
%                     B=[100;-40;0];
%                     C=[100;-40;20];
%                     D=[100;0;20];
                    P = [B,A;C,D];
                    X = P([1,4],:);
                    Y = P([2,5],:);
                    Z = P([3,6],:);
                    hh = surf(X,Y,Z);
                    hold on;
                    set(hh,'FaceColor','b');
                    
                    %%%%%%%%%-----�ж��ɻ����ز��������Ƿ�����ཻ����������������
                    if pp      %���ѡ����зɻ�λ�÷���
                        distan=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);
                        dirVector=[fly_x-CSB_X,fly_y-CSB_Y,fly_z-CSB_Z]/distan;

                        
                        
%                        syms x y z p %����3��������ע���ÿո�ָ����Һ��治�÷ֺţ�
%                        
%                        line_eqn1=(x-CSB_X)/(fly_x-CSB_X)==(y-CSB_Y)/(fly_y-CSB_Y);
%                        line_eqn2= (y-CSB_Y)/(fly_y-CSB_Y)==(z-CSB_Z)/(fly_z-CSB_Z);
%                        line_eqn=[x y z]==[CSB_X,CSB_Y,CSB_Z]+p.*dirVector; %�µ�����
                       BA=B-A;
                       BC=B-C;
                       Nor_obs=cross(BA,BC);
%                        AA=Nor_obs(1);
%                        BB=Nor_obs(2);
%                        CC=Nor_obs(3);
%                        plane_eqn=AA*(x-B(1))+BB*(y-B(2))+CC*(z-B(3))==0;
%                        eqn=[line_eqn plane_eqn];
%                        var=[x y z p];
% %                         eqn = sin(x) == x^2 - 1;
%                        [solx, soly, solz,solt] = solve(eqn, var);
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
                       set(handles.edit_crosspoint,'String','���ϰ���ƽ�治�ཻ');
                       inplane=false;
                       
                   else
                       pp=((n1-m1)*vp1+(n2-m2)*vp2+(n3-m3)*vp3) / (vp1* v1+ vp2* v2+ vp3* v3);
                       x=m1+v1*pp;
                       y=m2+v2*pp;
                       z=m3+v3*pp;
                       plot3(x,y,z,'o','Color',[1,0,0],'LineWidth',2); 
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
%                               msg=get(handles.edit_crosspoint,'String');
%                               str=[msg, 10];
%                               str=[str, 'ƽ�棺'];
%                               str=[str, num2str(i)];
%                               str=[str,  ',�ϰ���ƽ���ཻ�����ϰ�������'];
str1=['��������X��',num2str(x),10];
str2=['��������Y��',num2str(y),10];
str3=['��������Z��',num2str(z),10];
str=[str1,str2,str3];

                              set(handles.edit_crosspoint,'String',str);
                          else
                              inplane=false;
                               set(handles.edit_crosspoint,'String','�ϰ���ƽ���ཻ���������ϰ�������');
                          end
                              
                          
                          
                       
                   end
                   
                       
                  
                    end
                    
                    
                    
                    
                    %%%%%%%%%%%------�ж��ɻ����ز��������Ƿ�����ཻ

end

end

% axis equal;
% set(gca,'xtick',V); ����X��Ŀ̶ȱ��λ�� ax.XTick = [2 4 6 8 10]
% set(gca,'xticklable',sss); ����X���Ƿ��Ż����� ax.XTickLabel = {'Jan','Feb','Mar','Apr'}
% Example: ax.XLim = [0 10]
% 
% Example: ax.YLim = [-inf 10]
% 
% Example: ax.ZLim = [0 inf]
% 
% Alternatively, use the xlim, ylim, and zlim functions to set the limits.
% ax=axes(handles.axes1);    %ȡ����ͼ������
xmax=get(gca,'XLim');
ymax=get(gca,'YLim');
plot3([0,0],[0,ymax(2)]*0.8,[0,0],'LineWidth',2);
quiver3(0,ymax(2)*0.8,0,0,0.05*ymax(2),0,'LineWidth',2);%[0,0],[0,0],[0 0],[0 0],[0 ymax(2)*0.9],[0 0],'LineWidth',2);     %��һ����ͷ��Y��
quiver3([0,0],[0,0],[0,0],[0,0],[0,0],[0,5*h]);                    %��һ����ͷ��Z��
quiver3([0,0],[0,0],[0 0],[0 xmax(2)*0.9],[0 0],[0 0],'LineWidth',2);     %��һ����ͷ��X��
view(-11,34);   % ���ú��ʵ��ӽǡ�
% annotation('arrow',[0,0.7],[0,0.8]);


Ptable=PointTable;




% p = [0 0 0;
%     1 1 1;
%     3 2 1];
% r = diff(p);
% plot3(p(:,1),p(:,2),p(:,3),'+')
% hold on
% syms x y z
% % fsurf(solve(cross(r(1,:),r(2,:))*([x;y;z]-p(1,:)'),z))
% hold off

% plot3(0,0,-h,'*')

% plot(t,sin(2*t),'-mo',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)


%  axis([-10*r1 1000 -10*r1 10*r1 -10*r1 10*r1]);

% 
% h=0;  %�߶�
% r=handles.DVOR{2};
% pos=[0,0];
% t=0:0.0001:(2*pi);
% t=[t,0];
% xx=pos(1)+r*sin(t);
% yy=pos(2)+r*cos(t);
% zz=h*ones(size(t));
% % plot3(xx,yy,zz);
% mesh(xx,yy,zz);
% axis([-1.5*r 1.5*r -1.5*r 1.5*r -1.5*r 1.5*r]);
% % axis square;
% % set(gca,'XTick',0:10:100);
% % set(gca,'YTick',0:10:100);







% 
% t = 0:pi/20:2*pi;
% plot(t,sin(2*t),'-mo',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)
% Z = magic(5);
% b = bar3(Z);
% colorbar
% for k = 1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'interp';
% end

%%
%%




% --- Outputs from this function are returned to the command line.
function varargout = DVOR3DModel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

draw_DVOR_BASIC(hObject,eventdata,handles);
% global figOption;
%  global plotCount;
%  
% axes(handles.axes1);
% 
% 
% 
% 
% 
% 
% 
% % % % 
% % % % %%%%%%%%%%%%%%%%��ȡ�ϰ�������������ϰ���ƽ�棬���λ򳤷��Σ������ϱ߿��±߿��߶ȣ���б����ת�Ȳ�������
% % % %    dt1{1,1}=500;   %����
% % % %                         dt1{1,2}=30;    %��λ
% % % %                         dt1{1,3}=40;    %�ϱ߿�
% % % %                         dt1{1,4}=60;    %�±߶�
% % % %                         dt1{1,5}=80;    %�߶ȣ���ԭʼ�߶ȵ��ڸ�ֵ����TAN(��ǣ�
% % % %                       
% % % %                         dt1{1,6}=5;     %%Z������ƫ�� %3.3f
% % % %                         dt1{1,7}=20;    %Z����ת
% % % %                         dt1{1,8}=3;     %��б�Ƕ�
% % % %                         dt1{1,9}=1.0;   %����ϵ��
% % % %                         dt1{1,10}=true; %���ñ��
% 
% dt1=get(handles.uitable1,'Data');
% [a,b]=size(dt1);
% for i=1:a
%                        dist=dt1{i,1};   %����,��ͶӰ��XOYƽ��ľ��룬б�����
%                        az= dt1{i,2};    %��λ
%                         top_width=dt1{i,3};    %�ϱ߿�
%                         button_width=dt1{i,4};    %�±߶�
%                         elev_h=dt1{i,5};    %�߶ȣ���ԭʼ�߶ȵ��ڸ�ֵ����TAN(��ǣ�
%                       
%                         Z_shift=dt1{i,6};     %%Z������ƫ��,��
%                         Z_rotage=dt1{i,7};    %Z����ת
%                         Z_climing=dt1{i,8};     %��б�Ƕ�
%                         reflect_factor=dt1{i,9};   %����ϵ��
%                         active=dt1{i,10}; %���ñ��
%  
% if active==true     %%�Ƚ�ƽ���Ƶ�������ԭ�㴦��������x=0��ƽ���ϣ�������Ϊy����,��Ϊx����,ֽ�ⷽ��ΪZ����
%                     origin_H=elev_h/cos(Z_climing*pi/180);
%                     A0=[0;button_width/2;0];
%                     B0=[0;-button_width/2;0];
%                     C0=[0;-top_width/2;origin_H];
%                     D0=[0;top_width/2;origin_H];
%                      
%                     shift_T=[dist*sin(az*pi/180);dist*cos(az*pi/180);Z_shift];
%                     shift_y1=[0;button_width/2;0];
%                     shift_y2=[0;top_width/2;0];
%                     A=A0+shift_T;
%                     B=B0+shift_T;
%                     C=C0+shift_T;
%                     D=D0+shift_T;
%                     
%     
% %                     A=[100;0;0];
% %                     B=[100;-40;0];
% %                     C=[100;-40;20];
% %                     D=[100;0;20];
%                     P = [B,A;C,D];
%                     X = P([1,4],:);
%                     Y = P([2,5],:);
%                     Z = P([3,6],:);
%                     h = surf(X,Y,Z);
%                     hold on;
%                     set(h,'FaceColor','b');
% 
% end
% 
% end


% p = [0 0 0;
%     1 1 1;
%     3 2 1];
% r = diff(p);
% plot3(p(:,1),p(:,2),p(:,3),'+')
% hold on
% syms x y z
% % fsurf(solve(cross(r(1,:),r(2,:))*([x;y;z]-p(1,:)'),z))
% hold off

% plot3(0,0,-h,'*')

% plot(t,sin(2*t),'-mo',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)


%  axis([-10*r1 1000 -10*r1 10*r1 -10*r1 10*r1]);

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually havea white background on Windows.
%       See ISPC and COMP 



if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OpenFiles.
function pushbutton_OpenFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global delimiter;
global txtencoding;
% if txtencoding==''
%     txtencoding='UTF-8';
% end

[filename,filepath,fileindex]=uigetfile({'*.txt','�ı��ļ�(*.txt)';'*.xls','Excel�ļ�(*.xls)';'*.xlsx','Excel�ļ�(*.xlsx)';'*.csv','CSV�ļ�(*.csv)';'*.*','�����ļ�(*.*)'},'ѡ���ļ�');  %���ļ��Ի���

            strfile=[filepath filename];   %ƴ�Ӿ���·��
            set(handles.edit1,'String',strfile);  %

if fileindex~=0   %���û�е��ȡ��
    L1=length(filename);   %filename������׺��
    if L1<=4
        errordlg('�����ļ�','�ļ��򿪴���');  %����Ի���
        return;
    end
%     test=filename(1,L1-3:L1);   %�ļ�������ȡ��׺��
%     i = find('.'==filename);
% %ȥ���ļ���׺����ȡ�������ļ���
% imname = filename(1:i-1);
   [pathstr,name,test]=fileparts(strfile); % pathstr ���Ϊ E:\test ;  name ���Ϊ  test;suffix ���Ϊ  .txt

   set(handles.edit_figureTitle,'String',name);  %
   
    switch test   % test�����򿪵��ļ�����
        case  {'.xls','.XLS'}
            h=waitbar(0,'���ڶ�ȡ�ļ�....');   %������(���ļ��Ƚ���)
            [dataArray,txtdata,rawdata]=xlsread(strfile);  %���ݾ���·������xls�ļ���chengji��һ����ά��������xingming��һ��һά����������һ������chengjiֻ���洿���֣��ڶ�������xingmingֻ�����ַ�����
            waitbar(1,h,'���');
            delete(h);  
            
            [rrr,ccc]=size(rawdata);
             delimiter_len=ccc;
      set(handles.XCol_List,'String',1:delimiter_len);
       set(handles.YCol_List,'String',1:delimiter_len);
       set(handles.ZCol_List,'String',1:delimiter_len);
            set(handles.YCol_List,'Value',2);
            
            outtable= rawdata;
            set(handles.uitable1,'Data',outtable);
            handles.dTable=outtable;
           
           
            guidata(hObject, handles);   %����gui���� 
            
%            case  {'.csv','.CSV'}     %����Ǳ���ļ�
% 
%             h=waitbar(0,'���ڶ�ȡ�ļ�....');   %������(���ļ��Ƚ���)
%            rawdata=csvread(strfile);
% %             [chengji, xingming]=xlsread(strfile);  %���ݾ���·������xls�ļ���chengji��һ����ά��������xingming��һ��һά����������һ������chengjiֻ���洿���֣��ڶ�������xingmingֻ�����ַ�����
%             waitbar(1,h,'���');
%             delete(h);    
%              
%             [rrr,ccc]=size(rawdata);
%              delimiter_len=ccc;
%       set(handles.XCol_List,'String',1:delimiter_len);
%        set(handles.YCol_List,'String',1:delimiter_len);
%        set(handles.ZCol_List,'String',1:delimiter_len);
%        
%             outtable= rawdata;
%             set(handles.uitable1,'Data',outtable);
%             handles.dTable=outtable;
%            
%            
%             guidata(hObject, handles);   %����gui���� 
            
        case  {'.XLSX','.xlsx'}    %����Ǳ���ļ�
%             str=[filepath filename];   %ƴ�Ӿ���·��
%             set(handles.edit1,'String',str);  %
            h=waitbar(0,'���ڶ�ȡ�ļ�....');   %������(���ļ��Ƚ���)
           [dataArray,txtdata,rawdata]=xlsread(strfile);  %���ݾ���·������xls�ļ���chengji��һ����ά��������xingming��һ��һά����������һ������chengjiֻ���洿���֣��ڶ�������xingmingֻ�����ַ�����
            waitbar(1,h,'���');
            delete(h);    
%              data=[dataArray{1:end-1}];
            
            [rrr,ccc]=size(rawdata);
             delimiter_len=ccc;
      set(handles.XCol_List,'String',1:delimiter_len);
       set(handles.YCol_List,'String',1:delimiter_len);
       set(handles.ZCol_List,'String',1:delimiter_len);
         set(handles.YCol_List,'Value',2);
         
            outtable=rawdata;
            set(handles.uitable1,'Data',outtable);
            handles.dTable=outtable;
           
           
            guidata(hObject, handles);   %����gui����  
            
        case  {'.txt','.TXT','.csv','.CSV'}
%             strfile=[filepath filename];   %ƴ�Ӿ���·��
%             set(handles.edit1,'String',strfile);  
              h=waitbar(0,'���ڶ�ȡ�ļ�....');   %������(���ļ��Ƚ���)
%             fin=fopen(strfile,'r');   %��txt�ļ�
%             str=fgetl(fin);     %��ȡtxt�ļ��ĵ�һ��
            %   A=fscanf(fin,'%d','HeaderLines',1);  ����ǰһ�еı�������  
            %   [A B]=textscanf(fin,'%d %d');   ����һ��cell���͵�����[A B]��һ��һ�еĶ�  
            %   ftell();    �õ���ȡ��λ��
            
%             [str1 str2 str3 str4]=strread(str,'%s %s %s %s','delimiter',' '); %�Կո�ָ��ȡÿһ��
%             xingming(1)=str1;
%             counter=2;  %txt�ļ��ĵ�һ����(name yuwen shuxue yuwen)�����Դӵڶ��в�����Ҫ�����ݡ�
%             while feof(fin)==0   %����ܶ���txt�ļ������ݣ�(û�е�txt�ļ��Ľ�β)
%                 str=fgetl(fin);  %��ȡtxt�ļ���һ������
%                 [name yuwen shuxue yingyu]=strread(str,'%s %d %d %d','delimiter',' ');
%                 xingming(counter)=name;
%                 chengji(counter-1,:)=[yuwen shuxue yingyu];
%                 counter=counter+1;
%             end








%% ����������Ϊ�ı���ȡ:
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
%formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';  %'47������ 
formatSpec='';

%% ���ı��ļ���
fileID = fopen(strfile,'r','n',txtencoding);  %���ݱ����ʽ��ȡTXT�ļ�
           str=fgetl(fileID);     %��ȡtxt�ļ��ĵ�һ��
%            [filename, permission, machineformat, encoding] = fopen(fileID);

     if strcmp(test,'.csv') || strcmp(test,'.CSV')  % �ж��Ƿ�CSV�ļ���CSV�ļ���TXT�ļ��ķָ�����һ����CSV�Ƕ��ţ�
         delimiter=',';
     end
     
      CC=strsplit(str,delimiter);  %�ɵ�һ�У�ͨ���ָ������õ��ı��ļ���������
      
      delimiter_len=length(CC);
      set(handles.XCol_List,'String',1:delimiter_len);
       set(handles.YCol_List,'String',1:delimiter_len);
       set(handles.ZCol_List,'String',1:delimiter_len);
       
       set(handles.YCol_List,'Value',2);
       
      for          i=1:1:delimiter_len
        formatSpec=strcat(formatSpec,'%s') ;
      end
        formatSpec=strcat(formatSpec,'%[^\n\r]') ;  
%% ���ݸ�ʽ��ȡ�����С�
% �õ��û������ɴ˴������õ��ļ��Ľṹ����������ļ����ִ����볢��ͨ�����빤���������ɴ��롣
% % %   frewind(fileID)
% % % This call is identical to
% % % 
% % % fseek(fileID, 0, 'bof')
% % % Extended Capabilities
 frewind(fileID); %���ļ�ָָ�븴λ
 %%
 % 
 % <<FILENAME.PNG>>
 % 
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% �ر��ı��ļ���
fclose(fileID);
% hline=get(handles.checkbox_headline,'Value');
% if hline
% data=[dataArray{2:end-1}];
% else
    data=[dataArray{1:end-1}];
% end
  waitbar(1,h,'���');
     delete(h);    %���������
% raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
%% ��������ֵ�ı���������ת��Ϊ��ֵ��
% ������ֵ�ı��滻Ϊ NaN��
% % % % % % % raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% % % % % % % for col=1:length(dataArray)-1
% % % % % % %     raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
% % % % % % % end
% % % % % % % numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% % % % % % % 
% % % % % % % for col=1:1:delimiter_len   %[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
% % % % % % %     % ������Ԫ�������е��ı�ת��Ϊ��ֵ���ѽ�����ֵ�ı��滻Ϊ NaN��
% % % % % % %     rawData = dataArray{col};
% % % % % % %     for row=1:size(rawData, 1)
% % % % % % %         % ����������ʽ�Լ�Ⲣɾ������ֵǰ׺�ͺ�׺��
% % % % % % %         regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
% % % % % % %         try
% % % % % % %             result = regexp(rawData(row), regexstr, 'names');
% % % % % % %             numbers = result.numbers;
% % % % % % %             
% % % % % % %             % �ڷ�ǧλλ���м�⵽���š�
% % % % % % %             invalidThousandsSeparator = false;
% % % % % % %             if numbers.contains(',')
% % % % % % %                 thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
% % % % % % %                 if isempty(regexp(numbers, thousandsRegExp, 'once'))
% % % % % % %                     numbers = NaN;
% % % % % % %                     invalidThousandsSeparator = true;
% % % % % % %                 end
% % % % % % %             end
% % % % % % %             % ����ֵ�ı�ת��Ϊ��ֵ��
% % % % % % %             if ~invalidThousandsSeparator
% % % % % % %                 numbers = textscan(char(strrep(numbers, ',', '')), '%f');
% % % % % % %                 numericData(row, col) = numbers{1};
% % % % % % %                 raw{row, col} = numbers{1};
% % % % % % %             end
% % % % % % %         catch
% % % % % % %             raw{row, col} = rawData{row};
% % % % % % %         end
% % % % % % %     end
% % % % % % % end
% % % % % % % 
% % % % % % % 
% % % % % % % %% ������ֵԪ���滻Ϊ NaN
% % % % % % % R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % ���ҷ���ֵԪ��
% % % % % % % raw(R) = {NaN}; % �滻����ֵԪ��

%% �����������
% 'outtable = cell2mat(raw);
outtable=cellstr(data);


%% �����ʱ����
% clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;



            set(handles.uitable1,'Data',outtable);
            handles.dTable=outtable;
            guidata(hObject, handles);
            %fclose(fin);  %�ر��ļ���
            
        
        clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

        otherwise
            errordlg('�ļ����ʹ���','�ļ�����');  
            return;
    end
%--------------------- 
end


% --- Executes on selection change in XCol_List.
function XCol_List_Callback(hObject, eventdata, handles)
% hObject    handle to XCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XCol_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XCol_List


% --- Executes during object creation, after setting all properties.
function XCol_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in YCol_List.
function YCol_List_Callback(hObject, eventdata, handles)
% hObject    handle to YCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns YCol_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from YCol_List


% --- Executes during object creation, after setting all properties.
function YCol_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_XLable_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XLable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XLable as text
%        str2double(get(hObject,'String')) returns contents of edit_XLable as a double


% --- Executes during object creation, after setting all properties.
function edit_XLable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XLable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_YLable_Callback(hObject, eventdata, handles)
% hObject    handle to edit_YLable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_YLable as text
%        str2double(get(hObject,'String')) returns contents of edit_YLable as a double


% --- Executes during object creation, after setting all properties.
function edit_YLable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_YLable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plotCount.
function popupmenu_plotCount_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plotCount contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plotCount


% --- Executes during object creation, after setting all properties.
function popupmenu_plotCount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plotCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_headline.
function checkbox_headline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_headline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_headline


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not createntil after all CreateFcns called
global txtencoding;
global delimiter;
global figOption;
global plotCount;
txtencoding='GBK';
delimiter='\t';
figOption=1;
plotCount=1;

function edit_figureTitle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_figureTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_figureTitle as text
%        str2double(get(hObject,'String')) returns contents of edit_figureTitle as a double


% --- Executes during object creation, after setting all properties.
function edit_figureTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_figureTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global delimiter;
switch get(hObject,'Tag')   
    case 'radiobutton_tab'
       delimiter='\t';
    case 'radiobutton_comma'
        delimiter=',';
    case 'radiobutton_space'
        delimiter=' ';
    case 'radiobutton_semicolon'
        delimiter=';';
        
end


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global txtencoding;

% 'Big5'
% 'ISO-8859-1'
% 'windows-847'
% 'Big5-HKSCS'
% 'ISO-8859-2'
% 'windows-949'
% 'CP949'
% 'ISO-8859-3'
% 'windows-1250'
% 'EUC-KR'
% 'ISO-8859-4'
% 'windows-1251'
% 'EUC-JP'
% 'ISO-8859-5'
% 'windows-1252'
% 'EUC-TW'
% 'ISO-8859-6'
% 'windows-1253'
% 'GB18030'
% 'ISO-8859-7'
% 'windows-1254'
% 'GB2312'
% 'ISO-8859-8'
% 'windows-1255'
% 'GBK'
% 'ISO-8859-9'
% 'windows-1256'
% 'IBM866'
% 'ISO-8859-11'
% 'windows-1257'
% 'KOI8-R'
% 'ISO-8859-13'
% 'windows-1258'
% 'KOI8-U'
% 'ISO-8859-15'
% 'US-ASCII'
%  	
% 'Macintosh'
% 'UTF-8'
%  	
% 'Shift_JIS'
% %  

switch get(hObject,'Tag')   
    case 'radiobutton_UTF8'
       txtencoding='UTF-8';
    case 'radiobutton_GBK'
        txtencoding='GBK';
    case 'radiobutton_ASCII'
        txtencoding='US-ASCII ';
    case 'radiobutton_Big5'
        txtencoding='Big5';
     case 'radiobutton_JIS'    
        txtencoding='Shift_JIS';
        case 'radiobutton_GB2312'    
        txtencoding='GB2312'; 
end


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%  hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
newData = get(hObject,'Data'); %��ȡ���ݾ���

hang = eventdata.Indices;  %��ȡ������
if ~isempty(hang)
hangIndex = hang(1);  %��������ֵ
handles.hangIndex = hangIndex;  %����������ӵ��ṹ��
guidata(hObject, handles);  %���½ṹ��

 





end
% --------------------- 


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to delEle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%arr=get(handles.uitable1,'Data')
hangIndex = handles.hangIndex;  %��ȡѡ���Ժ���� ������
newData = get(handles.uitable1,'Data');  %��ȡ������ݾ���
newData(hangIndex,:) = [];   %ɾ��ѡ�е�ĳ������
set(handles.uitable1,'Data',newData);  %��ʾ�������
 handles.dTable=newData;
 guidata(hObject, handles);  %���½ṹ��
% save('newData.mat','newData');  %ɾ���Ժ󣬱���һ������
% --------------------- 
 


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figOption;
switch get(hObject,'Tag')   
    case 'radiobutton_newFigure'
      figOption=1;
    case 'radiobutton_newPlot'
        figOption=2;
    case 'radiobutton_newCur'
       figOption=3;
    
end        


% --- Executes on selection change in ZCol_List.
function ZCol_List_Callback(hObject, eventdata, handles)
% hObject    handle to ZCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ZCol_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ZCol_List


% --- Executes during object creation, after setting all properties.
function ZCol_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZCol_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_CurFit.
function pushbutton_CurFit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CurFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.dTable;

hline=get(handles.checkbox_headline,'Value');
if hline
data=data(2:end,:);
 
end
 
XCol=get(handles.XCol_List,'Value');
YCol=get(handles.YCol_List,'Value');
% ZCol=get(handles.ZCol_List,'Value');

ftitle=get(handles.edit_figureTitle, 'String');  %ͼ�α�����

 

 xData=str2double(data(:,XCol));
yData=str2double(data(:,YCol));
% ZZ=str2double(data(:,ZCol));

[xData_1, yData_1] = prepareCurveData( xData, yData );

% Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Robust = 'Bisquare';

ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';


% Fit model to data.
[fitresult, gof] = fit( xData_1, yData_1, ft, opts );

% % Plot fit with data.
% figure( 'Name', ftitle );
% h = plot( fitresult, xData_1, yData_1 );
% legend( h, 'yData vs. xData', 'Cur fitting', 'Location', 'NorthEast' );
% % Label axes
% xlabel xData
% ylabel yData
% grid on

% legend( h, surfacename, 'y vs. x','Location', 'NorthEast' );
a1=strcat('б��p1=',num2str(fitresult.p1),'  �Ƕȣ�',num2str(atan(fitresult.p1)*180/pi()),'��');
a2=strcat('�ؾ�p2=',num2str(fitresult.p2));
axes(handles.axes1);
cla;
curname=strcat('Y=',formula(fitresult));
%  
% plot(XX,YY,'DisplayName',ftitle);
% % plot(x,y1,'DisplayName','sin(x)')
% legend('show');
%  grid on;
% grid minor;
% hold off;
 h1 = plot( fitresult, xData_1, yData_1);
 
%  legend( h1, surfacename, 'y vs. x','Location', 'NorthEast' );
 hold on;
  h2 = plot( 0,0 );
  hold on;
  h3=plot(0,0)
  
  lgd = legend('Y vs. X',curname, a1,a2);
title(lgd,strcat(ftitle,'    Cur Fitting Result'));

%   legend( h2, a1,a2,'Location', 'NorthEast' );
%   lgd = legend('show');
lgd.FontSize = 10;
lgd.TextColor = 'blue';
legend('boxoff');
  hold off;
%    h3 = plot( fitresult, xData_1, yData_1,'DisplayName',a1 );
%     h4 = plot( fitresult, xData_1, yData_1 ,'DisplayName',a2 );
%     
% 
%   strs={surfacename,;a1;a2 };
%   ax=axis();
%   
%  text(ax(2)/2,ax(4)*0.2,strs); %{'\bfS
 cftool(xData_1, yData_1);

% --- Executes on button press in pushbutton_SurfaceFit.
function pushbutton_SurfaceFit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SurfaceFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% data=handles.dTable;

hline=get(handles.checkbox_headline,'Value');
if hline
data=data(2:end,:);
 
end
 
XCol=get(handles.XCol_List,'Value');
YCol=get(handles.YCol_List,'Value');
ZCol=get(handles.ZCol_List,'Value');

ftitle=get(handles.edit_figureTitle, 'String');  %ͼ�α�����

 

XX=str2double(data(:,XCol));
YY=str2double(data(:,YCol));
ZZ=str2double(data(:,ZCol));
[xData, yData, zData] = prepareSurfaceData(XX, YY, ZZ );


% ft = 'biharmonicinterp';


% % Fit model to data.
% [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
% 
% % Plot fit with data.
% figure( 'Name', ftitle );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, ftitle, 'z vs. x, y', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% zlabel z
% grid on
% view( -8.4, 62.0 );


% Set up fittype and options.
ft = fittype( 'poly11' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

axes(handles.axes1);
cla;
surfacename=strcat('Z=',formula(fitresult));

a1=strcat('��׼��߳� p00=',num2str(fitresult.p00));
a2=strcat('�����¶�FSL p10=',num2str(fitresult.p10),' = ',num2str(atan(fitresult.p10)*180/pi()),'��');
a3=strcat('�����¶�SSL p01=',num2str(fitresult.p01),' =',num2str(atan(fitresult.p01)*180/pi()),'��');
% a1=strcat('б��p1=',num2str(fitresult.p1),'  ','�Ƕȣ�',num2str(atan(fitresult.p1)*180/pi()),'��');
% a2=strcat('�ؾ�p2=',num2str(fitresult.p2));
axes(handles.axes1);
cla;


 h = plot( fitresult, [xData, yData], zData );
 
 hold on;
  h2 = plot( 0 );
  hold on;
  h3=plot(0);
  hold on;
  h4=plot(0);
   hold on;
  h5=plot(0);
  
  
%   tttile={'Z vs. X��Y',surfacename, a1,a2,a3};
  lgd = legend('Z vs. X��Y',surfacename, a1,a2,a3);
title(lgd,strcat(ftitle,'    Surface Fitting Result'));
% title(tttile);

%   legend( h2, a1,a2,'Location', 'NorthEast' );
%   lgd = legend('show');
lgd.FontSize = 9;
lgd.TextColor = 'blue';
legend('boxoff');
  hold off;
 
 
% legend( h, surfacename, 'z vs. x, y','Location', 'NorthEast' );
% 
%   strs={surfacename,;a1;a2;a3};
%   
%    ax=axis();
%   
%  text(ax(2)/2,ax(4)*0.2,strs); 
 
 cftool(xData, yData, zData);
 
%  text(0,0,strs); %{'\bfSig=\alpha*{Pb}^\beta'; 'd=2  D=2 mm  R=0mm'};
 
% % Plot fit with data.
% figure( 'Name', ftitle;);
% h = plot( fitresult, [xData, yData], zData );
% legend( h, ftitle, 'z vs. x, y', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% zlabel z
% grid on
% view( 43.6, 22.5 );


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dt=get(handles.uitable1,'Data');
[a,b]=size(dt);
% dt2=[a,b];

for i=1:a
    for j=1:b
        dt2(i,j)=string(dt(i,j));
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
    end
end

      
%  dt2=cell2mat(dt);

% save uidata.mat dt2

[filename,filepath]=uiputfile('.txt','save file');
% tt=["���ߺ�","�Ƕ�(��)","����(��)","�߶�(��)","����(%)","��λ(��)","��Ч"];
if filename~=0
str=[filepath,filename];
fid=fopen(str,'wt');  % -tģʽ����ʾ���ı���ʽ�������Ƕ�����ģʽ���ж�д
% fprintf(fid,[repmat('%s\t',1,size(dt2,2)),'\n'],tt'); %
fprintf(fid,[repmat('%s\t',1,size(dt2,2)),'\n'],dt2'); %

fclose(fid);
% dlmwrite(str,dt2,'delimiter','\t','-append');
end


% --- Executes on button press in radiobutton_space.
function radiobutton_space_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_space


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(handles.ht);
clc;
clear;


% --- Executes during object creation, after setting all properties.
function uibuttongroup4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_Add.
function pushbutton_Add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt={'���루m)','��λ���㣩','�ϱ߿�m)','�±߿�m)','�߶�(m)','Z��ֱƫ��(m)','Z��ת(��)','��б(��)','����ϵ��'};
title='����������';
lines=[1,1,1,1,1,1,1,1,1];
def={'800','90','60','60','80','5','45','0','1.0'};
tab=inputdlg(prompt,title,lines,def);
r1=str2double(tab{1});
r2=str2double(tab{2});
r3=str2double(tab{3});
r4=str2double(tab{4});
r5=str2double(tab{5});
r6=str2double(tab{6});
r7=str2double(tab{7});
r8=str2double(tab{8});
r9=str2double(tab{9});

newArray={r1,r2,r3,r4,r5,r6,r7,r8,r9,true}; %cell2mat(tab);
% newArray=cat(1,newrow,true);
olddata=get(handles.uitable1,'Data');
newData=[olddata;newArray];
set(handles.uitable1,'Data',newData);



function edit_az_Callback(hObject, eventdata, handles)
% hObject    handle to edit_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_az as text
%        str2double(get(hObject,'String')) returns contents of edit_az as a double


% --- Executes during object creation, after setting all properties.
function edit_az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_h_Callback(hObject, eventdata, handles)
% hObject    handle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_h as text
%        str2double(get(hObject,'String')) returns contents of edit_h as a double


% --- Executes during object creation, after setting all properties.
function edit_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_distance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_distance as text
%        str2double(get(hObject,'String')) returns contents of edit_distance as a double


% --- Executes during object creation, after setting all properties.
function edit_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plane.
function checkbox_plane_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plane



function edit_crosspoint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crosspoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crosspoint as text
%        str2double(get(hObject,'String')) returns contents of edit_crosspoint as a double


% --- Executes during object creation, after setting all properties.
function edit_crosspoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crosspoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
