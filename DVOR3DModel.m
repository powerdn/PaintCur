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
 
handles.ht=timer;  %添加时钟；计时器，用于刷新画面
set(handles.ht,'ExecutionMode','FixedRate');
set(handles.ht,'Period',1);     % 每秒刷新一次
set(handles.ht,'TimerFcn',{@dispnow,handles});  %绑定计时器执行FUNCTION
% Update handles structure
guidata(hObject, handles);
start(handles.ht);           %启动计时器

handles.DVOR=varargin;  %输入的参数：  DVOR3DModel(CSB_H,antenna_D,counterpoint_R,couterpoint_H,sb_basic_data);  %传递载波天线方高，天线阵直径，地网半径，地网高度，边带物理参数）

 guidata(hObject, handles)
 %%%%初始化uiTable中的数据，以2个plane为例
 
  dt1=cell(2,7);
                        dt1{1,1}=500;   %距离
                        dt1{1,2}=30;    %方位
                        dt1{1,3}=40;    %上边宽
                        dt1{1,4}=60;    %下边度
                        dt1{1,5}=80;    %高度，则原始高度等于该值除以TAN(倾角）
                      
                        dt1{1,6}=5;     %%Z轴上下偏移 %3.3f
                        dt1{1,7}=20;    %Z轴旋转
                        dt1{1,8}=3;     %倾斜角度
                        dt1{1,9}=1.0;   %反射系数
                        dt1{1,10}=true; %启用标记
       
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
 
draw_DVOR_BASIC(hObject,eventdata,handles);  %画DVOR地网、天线位置

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


%%  画DVOR台
function draw_DVOR_BASIC(hObject,eventdata,handles)

 %生成边带天线坐标  1：天线号，2：角度，3：距离，4：高度，5：幅度，6：相位，7：启用
   %传递载波天线方高，天线阵直径，地网半径，地网高度，边带物理参数）

  %         1          2          3        4           5  (48*6)     
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% % %   DVOR3DModel(CSB_H,antenna_D,counterpoint_R,couterpoint_H,sb_basic_data,handles);  %传递载波天线方高，天线阵直径，地网半径，地网高度，边带物理参数,父窗口句柄）
  
if isempty( handles.DVOR)     
   
    return;
end
global figOption;
 global plotCount;
 global Ptable;
 
axes(handles.axes1);    %取得作图区域句柄
cla;                   %清除作业区


% n=[1 1 1];   %;?%法向量n
% r=handles.DVOR{2};;  %?%圆的半径为1
% c=[1 1 1];   %?%圆心的坐标
% theta=(0:2*pi/100:2*pi)'; %;?%theta角从0到2*pi
% a=cross(n,[1 0 0]); %?%n与i叉乘，求取a向量
% if ~any(a)  %?%如果a为零向量，将n与j叉乘
%  a=cross(n,[0 1 0]);
% end
% b=cross(n,a); %求取b向量
% a=a/norm(a); %单位化a向量
% b=b/norm(b); %单位化b向量
% 
% c1=c(1)*ones(size(theta,1),1);
% c2=c(2)*ones(size(theta,1),1);
% c3=c(3)*ones(size(theta,1),1);
% 
% x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
% y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
% z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%圆上各点的z坐标
% 
% plot3(x,y,z)
% xlabel('x轴')
% ylabel('y轴')
% zlabel('z轴')

h_parent=handles.DVOR{6};

N=200;                % 画圆的分门辨率， 360°/200
r1=handles.DVOR{3};   %反射网半径
h=handles.DVOR{4};    %反射网高度
sbants=handles.DVOR{5};  %边带天物理参数，角度，高度，径向距离等
sb=sbants(:,1:3);        %边带天线的坐标 x,y,z

rr=linspace(0,r1,N);   %  径向步进长度
theta=linspace(0,2*pi,N); %步进角度
[R,the]=meshgrid(rr,theta);  %生成距阵，相当于N个同心圆，被分成N个扇区，生成R距阵，和the距阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %        R=[0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
% %           [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%           [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%            [0 1/N 2/N.....1]*r1                the=[0 1/N 2/N...1]*2pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=0*R.^0;   %生成Z的0距阵
[x,y]=pol2cart(the,R);   %极坐标转直角坐标
mesh(x,y,Z,'edgecolor',[0.2,0.5,1],'EdgeAlpha',0.5);    % 画网络，即地网的圆
 %,'FaceColor',[1 1 1]);
hold on;
% plot3(0,0,-2*rr ,'.');
% plot3(0,0,2*rr,'.')
colormap Colorcube;

%colorMapNames={'Parula','Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter','Gray','Bone','Copper','Pink','Lines','Colorube',
%                 'Prism','flag','white'

r2=handles.DVOR{2}/2;   %天线阵直径，半径=D/2
rr=linspace(0,r2,N);    
theta=linspace(0,4*pi,N);  % 4pi,画2圈
[R,the]=meshgrid(rr,theta);
Z=0.101*R.^0;          %画天线阵的圆，Z浮在零平面上，不然颜色会重叠
[x,y]=pol2cart(the,R);   %极坐标转直角坐标
C=x.*y;
mesh(x,y,Z,'edgecolor','y','EdgeAlpha',0.1);  %,'FaceColor',[0,1,1]);
 colormap Spring;
 
 p1=[0,0,2*r1];
 p2=[0,0,-h];
 XX=[p1(1),p2(1)];
  YY=[p1(2),p2(2)];
   ZZ=[p1(3),p2(3)];
plot3(XX,YY,ZZ,...             %%%画一条通过中央载波天线，即圆心的垂直线，从Z轴2*r1到-ｈ，取到地面。DVOR的坐标原点定在反射地网的中心，即中央载波天线的基础点。
         'LineWidth',2,...
          'MarkerEdgeColor','r',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10);

[a,b]=size(sb);    %画48根边带天线，利用一个”0=“和一条垂直线表示。
for i=1:a
    plot3(sb(i,1),sb(i,2),sb(i,3),'o','Color',[0,1,0],'LineWidth',2);   %画一个“o'形状，代表ALFORD环形天线
    plot3([sb(i,1),sb(i,1)],[sb(i,2),sb(i,2)],[0,sb(i,3)],'LineWidth',1,'Color',[1 0 0]);  %画垂直线，代表边带天线杆
end


%%%%%%%%%%%%%画飞机位置----------------------------------
pp=get(handles.checkbox_plane,'Value');
CSB_X=0;
CSB_Y=0;
CSB_Z=sb(1,3);

if pp
  fly_z=str2double(get(handles.edit_h,'String'))-h;   %飞机高度减去地网高度，是Z轴高度
  FlySimulate_Circle=str2double(get(handles.edit_distance,'String'));
  ang=str2double(get(handles.edit_az,'String'));
  
    if fly_z<0
       
       set(handles.edit_cross,'String',"飞机在地网下方！请重新设置飞行模式和参数。");   %如果数值为负，则退出。
    return;
    end
    
    fly_d=sqrt(FlySimulate_Circle^2-fly_z^2);
%     fly_angle=atan(fly_z/fly_d)*180/pi;   %仿真飞机的仰角，单位是度°
    fly_x=fly_d*sin(ang*pi/180);
    fly_y=fly_d*cos(ang*pi/180);
   plot3(fly_x,fly_y,fly_z,'o','Color',[1,0,0],'LineWidth',10);  %画一符号表示飞机。
   plot3([0,fly_x],[0,fly_y],[CSB_Z,fly_z],'LineWidth',1,'Color',[1 0 0]);  %画载波天线与飞机的连线 
end
%%%%%%%%%%%%------------结束-画飞机位置


% % % 
% % % %%%%%%%%%%%%%%%%读取障碍物参数，生成障碍物平面，梯形或长方形，，由上边宽，下边宽，高度，倾斜，旋转等参数构成
% % %    dt1{1,1}=500;   %距离
% % %                         dt1{1,2}=30;    %方位
% % %                         dt1{1,3}=40;    %上边宽
% % %                         dt1{1,4}=60;    %下边度
% % %                         dt1{1,5}=80;    %高度，则原始高度等于该值除以TAN(倾角）
% % %                       
% % %                         dt1{1,6}=5;     %%Z轴上下偏移 %3.3f
% % %                         dt1{1,7}=20;    %Z轴旋转
% % %                         dt1{1,8}=3;     %倾斜角度
% % %                         dt1{1,9}=1.0;   %反射系数
% % %                         dt1{1,10}=true; %启用标记

dt1=get(handles.uitable1,'Data');  %读取障碍物表格中的数据表
[a,b]=size(dt1);
iii=1;
for i=1:a       %循环生成障碍物平面，每行一个，得到N组障碍物梯形平南的座标。并赋值于公共变量Ptable.
                       dist=dt1{i,1};   %距离,是投影在XOY平面的距离，斜距另计
                       az= dt1{i,2};    %方位
                        top_width=dt1{i,3};    %上边宽
                        bottom_width=dt1{i,4};    %下边度
                        elev_h=dt1{i,5};    %高度，则原始高度等于该值除以TAN(倾角）
                      
                        Z_shift=dt1{i,6};     %%Z轴上下偏移,米
                        Z_rotate=dt1{i,7};    %Z轴旋转
                        Z_inclination=dt1{i,8};     %Z方向的倾斜角度
                        reflect_factor=dt1{i,9};   %反射系数
                        active=dt1{i,10};          %启用标记
 
if active==true     %%先将平面移到坐标轴原点处，并处于x=0的平面上，坐标上为y正向,右为x正向,纸外方向为Z正向
                    origin_H=elev_h/cos(Z_inclination*pi/180);     %  计算平面不倾斜时的高度
                    A0=[0;bottom_width/2;0];     %创建距阵用的是分号，代表用列的行列式。 表示X Y Z .梯形（距形）的四个项点坐标，A点在上，B点在下，C点在B上方，D点在A上方
                    B0=[0;-bottom_width/2;0];
                    C0=[0;-top_width/2;origin_H];
                    D0=[0;top_width/2;origin_H];
                    
                    
                    %  Z轴方向倾斜，这时AB不变，CD往外、内倾斜，则A0、B0不变，实则绕Y轴旋转
                    shift_Y=[cos(Z_inclination*pi/180),0,sin(Z_inclination*pi/180);0,1,0;-sin(Z_inclination*pi/180),0,cos(Z_inclination*pi/180)];
                     
                    A0=shift_Y*A0;
                    B0=shift_Y*B0;
                    C0=shift_Y*C0;
                    D0=shift_Y*D0;
                    
                    
                    %绕Z轴旋转
                    shift_Z=[cos(Z_rotate*pi/180),-sin(Z_rotate*pi/180),0;sin(Z_rotate*pi/180),cos(Z_rotate*pi/180),0;0,0,1];
                    
                    A0=shift_Z*A0;
                    B0=shift_Z*B0;
                    C0=shift_Z*C0;
                    D0=shift_Z*D0;
                    
                    
                    
                    %平移到 dist,Az处
                    shift_T=[dist*sin(az*pi/180);dist*cos(az*pi/180);Z_shift];
                    shift_y1=[0;bottom_width/2;0];
                    shift_y2=[0;top_width/2;0];
                    
                    
                    
                    
                    
                    
                    A=A0+shift_T;
                    B=B0+shift_T;
                    C=C0+shift_T;
                    D=D0+shift_T;
                    
                    PointTable(iii,1)={[A,B,C,D]};   %添加到公共变量.这是元胞数组cell。包括四个顶点和反射系数。只有ACTIVE的障碍物才添加，不是ACTIVE的不添加
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
                    
                    %%%%%%%%%-----判读飞机与载波的连线是否距形相交，并给出交点坐标
                    if pp      %如果选择进行飞机位置分析
                        distan=sqrt((fly_x-CSB_X)^2+(fly_y-CSB_Y)^2+(fly_z-CSB_Z)^2);
                        dirVector=[fly_x-CSB_X,fly_y-CSB_Y,fly_z-CSB_Z]/distan;

                        
                        
%                        syms x y z p %定义3个变量，注意用空格分隔，且后面不用分号；
%                        
%                        line_eqn1=(x-CSB_X)/(fly_x-CSB_X)==(y-CSB_Y)/(fly_y-CSB_Y);
%                        line_eqn2= (y-CSB_Y)/(fly_y-CSB_Y)==(z-CSB_Z)/(fly_z-CSB_Z);
%                        line_eqn=[x y z]==[CSB_X,CSB_Y,CSB_Z]+p.*dirVector; %新点坐标
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
%  知直线L过点m（m1，m2，m3），且方向向量为VL（v1，v2，v3），
%平面P过点n（n1，n2，n3），且法线方向向量为VP（vp1，vp2，vp3），求得直线与平面的交点O的坐标（x，y，z）：
% % 将直线方程写成参数方程形式，即有：
% % % x = m1+ v1 * t
% % y = m2+ v2 * t (1)
% % z = m3+ v3 * t
% % 将平面方程写成点法式方程形式，即有：
% % vp1 * (x – n1) + vp2 * (y – n2) + vp3 * (z – n3) = 0 (2)
%t = ((n1-m1)*vp1+(n2-m2)*vp2+(n3-m3)*vp3) / (vp1* v1+ vp2* v2+ vp3* v3);
%如果（3）式中分母(vp1* v1+ vp2* v2+ vp3* v3)为0，则表示直线与平面平行，即直线与平面没有交点。
%求解出t后，然后将t代入式（1）即可求得交点O的坐标（x，y，z）。
                        
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
                       set(handles.edit_crosspoint,'String','与障碍物平面不相交');
                       inplane=false;
                       
                   else
                       pp=((n1-m1)*vp1+(n2-m2)*vp2+(n3-m3)*vp3) / (vp1* v1+ vp2* v2+ vp3* v3);
                       x=m1+v1*pp;
                       y=m2+v2*pp;
                       z=m3+v3*pp;
                       plot3(x,y,z,'o','Color',[1,0,0],'LineWidth',2); 
                       Pt=[x; y; z];
                       d_AB=norm(cross(B-A,Pt-A))/norm(B-A);   %计算交点到四条边的距离，这是公式。用叉乘得到三角形面积，再除以边，就得到点到直线的距离
                        d_BC=norm(cross(C-B,Pt-B))/norm(C-B); 
                         d_CD=norm(cross(D-C,Pt-C))/norm(D-C); 
                          d_DA=norm(cross(A-D,Pt-D))/norm(A-D); 
                          %如果4个三角形面积之和等于梯形的面积，则Pt点在梯形内，否则不在
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
%                               str=[str, '平面：'];
%                               str=[str, num2str(i)];
%                               str=[str,  ',障碍物平面相交且在障碍物面上'];
str1=['交点坐标X：',num2str(x),10];
str2=['交点坐标Y：',num2str(y),10];
str3=['交点坐标Z：',num2str(z),10];
str=[str1,str2,str3];

                              set(handles.edit_crosspoint,'String',str);
                          else
                              inplane=false;
                               set(handles.edit_crosspoint,'String','障碍物平面相交，但不在障碍物面上');
                          end
                              
                          
                          
                       
                   end
                   
                       
                  
                    end
                    
                    
                    
                    
                    %%%%%%%%%%%------判读飞机与载波的连线是否距形相交

end

end

% axis equal;
% set(gca,'xtick',V); 设置X轴的刻度标记位置 ax.XTick = [2 4 6 8 10]
% set(gca,'xticklable',sss); 设置X轴标记符号或文字 ax.XTickLabel = {'Jan','Feb','Mar','Apr'}
% Example: ax.XLim = [0 10]
% 
% Example: ax.YLim = [-inf 10]
% 
% Example: ax.ZLim = [0 inf]
% 
% Alternatively, use the xlim, ylim, and zlim functions to set the limits.
% ax=axes(handles.axes1);    %取得作图区域句柄
xmax=get(gca,'XLim');
ymax=get(gca,'YLim');
plot3([0,0],[0,ymax(2)]*0.8,[0,0],'LineWidth',2);
quiver3(0,ymax(2)*0.8,0,0,0.05*ymax(2),0,'LineWidth',2);%[0,0],[0,0],[0 0],[0 0],[0 ymax(2)*0.9],[0 0],'LineWidth',2);     %画一个箭头，Y轴
quiver3([0,0],[0,0],[0,0],[0,0],[0,0],[0,5*h]);                    %画一个箭头，Z轴
quiver3([0,0],[0,0],[0 0],[0 xmax(2)*0.9],[0 0],[0 0],'LineWidth',2);     %画一个箭头，X轴
view(-11,34);   % 设置合适的视角。
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
% h=0;  %高度
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
% % % % %%%%%%%%%%%%%%%%读取障碍物参数，生成障碍物平面，梯形或长方形，，由上边宽，下边宽，高度，倾斜，旋转等参数构成
% % % %    dt1{1,1}=500;   %距离
% % % %                         dt1{1,2}=30;    %方位
% % % %                         dt1{1,3}=40;    %上边宽
% % % %                         dt1{1,4}=60;    %下边度
% % % %                         dt1{1,5}=80;    %高度，则原始高度等于该值除以TAN(倾角）
% % % %                       
% % % %                         dt1{1,6}=5;     %%Z轴上下偏移 %3.3f
% % % %                         dt1{1,7}=20;    %Z轴旋转
% % % %                         dt1{1,8}=3;     %倾斜角度
% % % %                         dt1{1,9}=1.0;   %反射系数
% % % %                         dt1{1,10}=true; %启用标记
% 
% dt1=get(handles.uitable1,'Data');
% [a,b]=size(dt1);
% for i=1:a
%                        dist=dt1{i,1};   %距离,是投影在XOY平面的距离，斜距另计
%                        az= dt1{i,2};    %方位
%                         top_width=dt1{i,3};    %上边宽
%                         button_width=dt1{i,4};    %下边度
%                         elev_h=dt1{i,5};    %高度，则原始高度等于该值除以TAN(倾角）
%                       
%                         Z_shift=dt1{i,6};     %%Z轴上下偏移,米
%                         Z_rotage=dt1{i,7};    %Z轴旋转
%                         Z_climing=dt1{i,8};     %倾斜角度
%                         reflect_factor=dt1{i,9};   %反射系数
%                         active=dt1{i,10}; %启用标记
%  
% if active==true     %%先将平面移到坐标轴原点处，并处于x=0的平面上，坐标上为y正向,右为x正向,纸外方向为Z正向
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

[filename,filepath,fileindex]=uigetfile({'*.txt','文本文件(*.txt)';'*.xls','Excel文件(*.xls)';'*.xlsx','Excel文件(*.xlsx)';'*.csv','CSV文件(*.csv)';'*.*','所有文件(*.*)'},'选择文件');  %打开文件对话框

            strfile=[filepath filename];   %拼接绝对路径
            set(handles.edit1,'String',strfile);  %

if fileindex~=0   %如果没有点击取消
    L1=length(filename);   %filename包含后缀名
    if L1<=4
        errordlg('错误文件','文件打开错误');  %错误对话框
        return;
    end
%     test=filename(1,L1-3:L1);   %文件名，截取后缀名
%     i = find('.'==filename);
% %去除文件后缀，提取单纯的文件名
% imname = filename(1:i-1);
   [pathstr,name,test]=fileparts(strfile); % pathstr 结果为 E:\test ;  name 结果为  test;suffix 结果为  .txt

   set(handles.edit_figureTitle,'String',name);  %
   
    switch test   % test是所打开的文件类型
        case  {'.xls','.XLS'}
            h=waitbar(0,'正在读取文件....');   %进度条(打开文件比较慢)
            [dataArray,txtdata,rawdata]=xlsread(strfile);  %根据绝对路径，打开xls文件。chengji是一个三维列向量，xingming是一个一维列向量。第一个参数chengji只保存纯数字，第二个参数xingming只保存字符串。
            waitbar(1,h,'完成');
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
           
           
            guidata(hObject, handles);   %更新gui数据 
            
%            case  {'.csv','.CSV'}     %如果是表格文件
% 
%             h=waitbar(0,'正在读取文件....');   %进度条(打开文件比较慢)
%            rawdata=csvread(strfile);
% %             [chengji, xingming]=xlsread(strfile);  %根据绝对路径，打开xls文件。chengji是一个三维列向量，xingming是一个一维列向量。第一个参数chengji只保存纯数字，第二个参数xingming只保存字符串。
%             waitbar(1,h,'完成');
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
%             guidata(hObject, handles);   %更新gui数据 
            
        case  {'.XLSX','.xlsx'}    %如果是表格文件
%             str=[filepath filename];   %拼接绝对路径
%             set(handles.edit1,'String',str);  %
            h=waitbar(0,'正在读取文件....');   %进度条(打开文件比较慢)
           [dataArray,txtdata,rawdata]=xlsread(strfile);  %根据绝对路径，打开xls文件。chengji是一个三维列向量，xingming是一个一维列向量。第一个参数chengji只保存纯数字，第二个参数xingming只保存字符串。
            waitbar(1,h,'完成');
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
           
           
            guidata(hObject, handles);   %更新gui数据  
            
        case  {'.txt','.TXT','.csv','.CSV'}
%             strfile=[filepath filename];   %拼接绝对路径
%             set(handles.edit1,'String',strfile);  
              h=waitbar(0,'正在读取文件....');   %进度条(打开文件比较慢)
%             fin=fopen(strfile,'r');   %打开txt文件
%             str=fgetl(fin);     %读取txt文件的第一行
            %   A=fscanf(fin,'%d','HeaderLines',1);  忽略前一行的标题内容  
            %   [A B]=textscanf(fin,'%d %d');   返回一个cell类型的数据[A B]，一列一列的读  
            %   ftell();    得到读取的位置
            
%             [str1 str2 str3 str4]=strread(str,'%s %s %s %s','delimiter',' '); %以空格分割，截取每一行
%             xingming(1)=str1;
%             counter=2;  %txt文件的第一行是(name yuwen shuxue yuwen)，所以从第二行才是需要的内容。
%             while feof(fin)==0   %如果能读到txt文件的内容，(没有到txt文件的结尾)
%                 str=fgetl(fin);  %读取txt文件的一行内容
%                 [name yuwen shuxue yingyu]=strread(str,'%s %d %d %d','delimiter',' ');
%                 xingming(counter)=name;
%                 chengji(counter-1,:)=[yuwen shuxue yingyu];
%                 counter=counter+1;
%             end








%% 将数据列作为文本读取:
% 有关详细信息，请参阅 TEXTSCAN 文档。
%formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';  %'47列数据 
formatSpec='';

%% 打开文本文件。
fileID = fopen(strfile,'r','n',txtencoding);  %根据编码格式读取TXT文件
           str=fgetl(fileID);     %读取txt文件的第一行
%            [filename, permission, machineformat, encoding] = fopen(fileID);

     if strcmp(test,'.csv') || strcmp(test,'.CSV')  % 判读是否CSV文件，CSV文件与TXT文件的分隔符不一样，CSV是逗号，
         delimiter=',';
     end
     
      CC=strsplit(str,delimiter);  %由第一行，通过分隔符，得到文本文件的列数。
      
      delimiter_len=length(CC);
      set(handles.XCol_List,'String',1:delimiter_len);
       set(handles.YCol_List,'String',1:delimiter_len);
       set(handles.ZCol_List,'String',1:delimiter_len);
       
       set(handles.YCol_List,'Value',2);
       
      for          i=1:1:delimiter_len
        formatSpec=strcat(formatSpec,'%s') ;
      end
        formatSpec=strcat(formatSpec,'%[^\n\r]') ;  
%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
% % %   frewind(fileID)
% % % This call is identical to
% % % 
% % % fseek(fileID, 0, 'bof')
% % % Extended Capabilities
 frewind(fileID); %将文件指指针复位
 %%
 % 
 % <<FILENAME.PNG>>
 % 
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% 关闭文本文件。
fclose(fileID);
% hline=get(handles.checkbox_headline,'Value');
% if hline
% data=[dataArray{2:end-1}];
% else
    data=[dataArray{1:end-1}];
% end
  waitbar(1,h,'完成');
     delete(h);    %清除进度条
% raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
%% 将包含数值文本的列内容转换为数值。
% 将非数值文本替换为 NaN。
% % % % % % % raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
% % % % % % % for col=1:length(dataArray)-1
% % % % % % %     raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
% % % % % % % end
% % % % % % % numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% % % % % % % 
% % % % % % % for col=1:1:delimiter_len   %[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
% % % % % % %     % 将输入元胞数组中的文本转换为数值。已将非数值文本替换为 NaN。
% % % % % % %     rawData = dataArray{col};
% % % % % % %     for row=1:size(rawData, 1)
% % % % % % %         % 创建正则表达式以检测并删除非数值前缀和后缀。
% % % % % % %         regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
% % % % % % %         try
% % % % % % %             result = regexp(rawData(row), regexstr, 'names');
% % % % % % %             numbers = result.numbers;
% % % % % % %             
% % % % % % %             % 在非千位位置中检测到逗号。
% % % % % % %             invalidThousandsSeparator = false;
% % % % % % %             if numbers.contains(',')
% % % % % % %                 thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
% % % % % % %                 if isempty(regexp(numbers, thousandsRegExp, 'once'))
% % % % % % %                     numbers = NaN;
% % % % % % %                     invalidThousandsSeparator = true;
% % % % % % %                 end
% % % % % % %             end
% % % % % % %             % 将数值文本转换为数值。
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
% % % % % % % %% 将非数值元胞替换为 NaN
% % % % % % % R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 查找非数值元胞
% % % % % % % raw(R) = {NaN}; % 替换非数值元胞

%% 创建输出变量
% 'outtable = cell2mat(raw);
outtable=cellstr(data);


%% 清除临时变量
% clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;



            set(handles.uitable1,'Data',outtable);
            handles.dTable=outtable;
            guidata(hObject, handles);
            %fclose(fin);  %关闭文件流
            
        
        clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

        otherwise
            errordlg('文件类型错误','文件错误');  
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
newData = get(hObject,'Data'); %获取数据矩阵

hang = eventdata.Indices;  %获取行索引
if ~isempty(hang)
hangIndex = hang(1);  %行索引赋值
handles.hangIndex = hangIndex;  %把行索引添加到结构体
guidata(hObject, handles);  %更新结构体

 





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
hangIndex = handles.hangIndex;  %获取选择以后传入的 行索引
newData = get(handles.uitable1,'Data');  %获取表格数据矩阵
newData(hangIndex,:) = [];   %删除选中的某行数据
set(handles.uitable1,'Data',newData);  %显示到表格中
 handles.dTable=newData;
 guidata(hObject, handles);  %更新结构体
% save('newData.mat','newData');  %删除以后，保存一次数据
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

ftitle=get(handles.edit_figureTitle, 'String');  %图形标题名

 

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
a1=strcat('斜率p1=',num2str(fitresult.p1),'  角度：',num2str(atan(fitresult.p1)*180/pi()),'°');
a2=strcat('截距p2=',num2str(fitresult.p2));
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

ftitle=get(handles.edit_figureTitle, 'String');  %图形标题名

 

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

a1=strcat('基准点高程 p00=',num2str(fitresult.p00));
a2=strcat('纵向坡度FSL p10=',num2str(fitresult.p10),' = ',num2str(atan(fitresult.p10)*180/pi()),'°');
a3=strcat('横向坡度SSL p01=',num2str(fitresult.p01),' =',num2str(atan(fitresult.p01)*180/pi()),'°');
% a1=strcat('斜率p1=',num2str(fitresult.p1),'  ','角度：',num2str(atan(fitresult.p1)*180/pi()),'°');
% a2=strcat('截距p2=',num2str(fitresult.p2));
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
  
  
%   tttile={'Z vs. X，Y',surfacename, a1,a2,a3};
  lgd = legend('Z vs. X，Y',surfacename, a1,a2,a3);
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
% tt=["天线号","角度(°)","距离(米)","高度(米)","幅度(%)","相位(°)","有效"];
if filename~=0
str=[filepath,filename];
fid=fopen(str,'wt');  % -t模式，表示按文本形式，而不是二进制模式进行读写
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
prompt={'距离（m)','方位（°）','上边宽（m)','下边宽（m)','高度(m)','Z垂直偏移(m)','Z旋转(°)','倾斜(°)','反射系数'};
title='请输入数据';
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
