function raw_data=hmd_distortion_simple
%addpath(genpath('C:\Users\Yuta\Dropbox\work\Projects\MATLAB_TOOLS\bitbucket'));
addpath(genpath('subfunctions'));

clear all
rng('default');
rng(512034961);

%return
USE_NEW_DATAFORMAT=true;
VISUALIZE_RAW_DATA=false;


if USE_NEW_DATAFORMAT==false
    file_K = 'C:\Users\Yuta\Dropbox\work\Projects\20140530_EyeHMDCalibDataset\calib_ids\cameraintrinsics.txt';
    root = 'C:\Users\Yuta\Dropbox\work\Projects\20140530_EyeHMDCalibDataset\calib_hmd_optics\' % Training data set path
    N_id = 23;% The maximum dataset index
    idset=[0:N_id];
    idset=[1:3 5:7 9:13 20:23]%
    %idset=[1:3 5:7 9:13 14:19 20:23]%
    idset=[1:3 5:7 9:13 20:23 24:27]%
    
    w_hmd=1280;
    h_hmd=1024;
    w_userview=1280;
    h_userview=1024;
else
    data_set = 1
    switch data_set
        case 1
            %kRoot = 'C:\Users\Yuta\Copy\20141214_EyeHMDCalibDataset_ST60_PS\';
            %kRoot='C:\Users\Yuta\Documents\ProjectLocal\20141214_EyeHMDCalibDataset_ST60_PS\';
            kRoot='.\';
            file_K = [kRoot,'data\cameraintrinsics.txt'];
            file_screen = [kRoot,'data\virtual_screen'];
            root   = [kRoot,'data\']; % Training data set path
            N_id = 12;% The maximum dataset index
            idset=[1:4 6:N_id]%
            N_id = 17
            TerminateNum=18;
            %idset=[18 1:N_id]%
            idset=[1:N_id]%
%             % VR demo 1
%             TerminateNum=19;
%             idset=[19]%
%             % VR demo 2
%             TerminateNum=20;
%             idset=[20]%
            w_hmd=1280;
            h_hmd=1024;
            w_userview=1280;
            h_userview=1024;
        case 2 % BT100, ISMAR2015
            kRoot = 'C:\Users\Yuta\Documents\ProjectLocal\20150123_HMDCalibDataset_BT100\';%C:\Users\Yuta\Dropbox\work\Projects\20150123_HMDCalibDataset_BT100\';
            file_K = [kRoot,'calib_ids\cameraintrinsics.txt'];
            root   = [kRoot,'DisplayCalibrationToolbox\data\']; % Training data set path
            N_id = 10;% The maximum dataset index
            idset=[1 3:6 9:N_id];
            %idset=[11];
            TerminateNum=11;       
            w_hmd=656;
            h_hmd=496;
            w_userview=2448;
            h_userview=2048;
    end
end

IGNORE_EYE_POSITION = true;
%IGNORE_EYE_POSITION = false;
NORMILIZE_DATA=true;
%NORMILIZE_DATA=false;
USE_IWLS = true;% Use importance-weighted regularized kernel least-squares
USE_IWLS = false;

DISTORT_UVST = true; % light field parameterization is UVST, i.e. lumigraph
%DISTORT_UVST = false;

% IDS intrinsic matrix
K_raw = openUbitrack3x3MatrixCalib(file_K);
%K=[1 1 0;0 1 0; 0 0 1];
%coordinates_type = analyse_intrinsic_matrix(K);
offsetw=w_hmd/2-0.5;
offseth=h_hmd/2-0.5;
offsetw=0;
offseth=0;
K_E = convert_intrinsic_matrix2OpenGL(K_raw,h_userview); %% convert in the OpenGL convention

subplot_row = 1;
subplot_col = 3;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load training datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%root = '20140525_optical_distortion_light_field\';
X=[]; % Observed 3D points
Y=[]; % Camera center
X_distorted=[]; % Distorted 3D points
XS_W=[]; % Points on the virtual plane in the world
XS_S=[]; % Points on the virtual plane in the screen coordinate system
XS_W_distorted=[]; % Points on the virtual plane in the world
XS_S_distorted=[]; % Points on the virtual plane in the screen coordinate system
UVST_distorted=[]; % UVST plane
UVST          =[]; % UVST plane
t_E0W=[];
G_E=[];
G_distorted_E=[];

%
misc.fig_id=1;
misc.xlim_min=400;
misc.xlim_max=1400;
misc.ylim_min=200;
misc.ylim_max=800;
misc.quiverscale=1;
misc.h = h_userview;
misc.USE_NEW_DATAFORMAT=USE_NEW_DATAFORMAT;

figure(misc.fig_id);clf;

% X0: 3D grid points in W
% Y0: 3D eye-camera position, t_E0W
    function [XS_W0,XS_S0] = intersetRayWithPlane(X0,Y0,t_PW,R_SW)
        R_WS = R_SW';
        t_WP = -R_SW'*t_PW;
        r3=R_WS(3,:)';
        dx=X0-Y0; % 3D grid points in E0
        a=( t_PW'*r3 - Y0'*r3)./(dx'*r3);% a: Nx1
        XS_W0= Y0+dx.*repmat(a',3,1);% 3D points on the virtual display in W
        XS_S0=R_WS*XS_W0+repmat(t_WP,1,size(XS_W0,2));% in S, so XS_S(3,:) are all 0
    end

%% VR paeper Fig
DRAW_VR_ONE=false;
if DRAW_VR_ONE
    figure(44);
    clf
    fig_width =900;
    fig_height=500;
    screen_size=get(gcf,'Position');
    set(gcf,'Position',[screen_size(1) screen_size(2) fig_width fig_height]);
    sub_fig_idx=2;
    sub_fig_num = length(idset);
    sub_fig_row = 4;
    sub_fig_col = ceil(sub_fig_num/sub_fig_row);
end
%%
figure(21);
clf;
figure(20);
clf;
for id =  idset
    if USE_NEW_DATAFORMAT
        sub = strcat(num2str(id),'\');
    else
        sub = strcat(num2str(id,'%.3d'),'\');
    end
    dir=strcat(root,sub);
    % X0: 3D grid points in W
    % Y0: 3D eye-camera position, t_E0W
    % X_distorted0: 3D grid points in W reprojected from 2D grid points in
    % the image coordinate system
    % G_E: 3D grid points in E
    % G_distorted_E: 3D distorted grid points in E
    [X0,Y0,X_distorted0,t_E0W0, t_W0E0, G_E0, G_distorted_E0]=loadDataSet2(dir,K_E, misc);
    
    %%% Copmute the intersection of the virtual screen and the 3D points in the eye coordinate system.
    %R_SW=R_WS';
    %t_SW=-R_WS'*t_WS;
    %t_SW_c=-R_WS'*t_WS_c;
    tmp=load(file_screen);    
    if 1 % use newly estimated virtual screen plane
        R_SW        = tmp.R_SW
        %t_SW        = tmp.t_SW_c;
        t_SW        = tmp.t_SW
        meter2pixel = tmp.alpha
        if 0 % CIC
        tmp=load('.\CIC_param\display_param_linear.mat');
        meter2pixel=tmp.alpha
        R_WS=tmp.R_WS;
        t_WS=tmp.t_WS;
        R_WS(2:3,:) = -R_WS(2:3,:);
        t_WS(2:3)   = -t_WS(2:3);
        R_SW=R_WS'
        t_SW=-R_WS'*t_WS
        t_SW=t_SW-R_SW*[(w_hmd/2-0.5);(h_hmd/2-0.5);0]/meter2pixel
        
        end
    end
    
    % find the intersectaions of light rays and the virtual display
    [XS_W0,XS_S0] = intersetRayWithPlane(X0,Y0,t_SW,R_SW);
    % for X_distorted0
    [XS_W_distorted0,XS_S_distorted0] = intersetRayWithPlane(X_distorted0,Y0,t_SW,R_SW);
    
    %% Light field in 4D uv-st representation
    t_SW_z0=t_SW;
    t_SW_z0(3)=0;
    % find the intersectaions of light rays and the virtual display
    [~,UV] = intersetRayWithPlane(X0,Y0,t_SW_z0,R_SW);
    % for X_distorted0
    [~,UV_distorted] = intersetRayWithPlane(X_distorted0,Y0,t_SW_z0,R_SW);
    UV = UV(1:2,:);
    UV_distorted = UV_distorted(1:2,:);
    ST = XS_S0(1:2,:);
    ST_distorted = XS_S_distorted0(1:2,:);
    UVST0 = [UV;ST];
    UVST_distorted0 = [UV_distorted;ST_distorted];
    % append data
    UVST           = [UVST UVST0];
    UVST_distorted = [UVST_distorted UVST_distorted0];
    if id==TerminateNum
        figure(1);
        clf;hold on;
        plotv(UVST0(3:4,:));
        plotv(UVST_distorted0(3:4,:));
        save('UVSTsets','UVST0','UVST_distorted0','X0','X_distorted0');
        return
    end
    
    %%% VR paeper Fig
    if DRAW_VR_ONE
        %%% Render input light field
        h0=figure(44);
        subplot(sub_fig_row,sub_fig_col,sub_fig_idx);
        box on;
        colormap('jet');
        draw_st_plane(UVST0,UVST_distorted0, meter2pixel,w_hmd,h_hmd)
        % remove tick labels of images inside the subplot
        if mod(sub_fig_idx,sub_fig_col)~=1 && sub_fig_idx~=2
            set(gca,'yticklabel',{[]})
        end
        if sub_fig_row~=ceil(sub_fig_idx/sub_fig_col)
            set(gca,'xticklabel',{[]})
        end
        sub_fig_idx=sub_fig_idx+1;
        sub_fig_col
    end
    %%%%
    
    if VISUALIZE_RAW_DATA % debug
        figure(20);
        cfg1 = 'b.';
        cfg2 = 'r.';
        subplot(subplot_row,subplot_col,1);
        hold on;grid on;% axis equal;
        plot(UV(1,:),ST(1,:),cfg1);
        plot(UV_distorted(1,:),ST_distorted(1,:),cfg2);
        legend('Ideal','Distorted');
        title('u-s plane');
        subplot(subplot_row,subplot_col,2);
        hold on;grid on; %axis equal;
        plot(UV(2,:),ST(2,:),cfg1);
        plot(UV_distorted(2,:),ST_distorted(2,:),cfg2);
        title('v-t plane');
        subplot(subplot_row,subplot_col,3);
        hold on;grid on; %axis equal;
        x0=ST(1,:)*meter2pixel;
        y0=ST(2,:)*meter2pixel;
        plot(x0,y0,cfg1);
        x1=ST_distorted(1,:)*meter2pixel;
        y1=ST_distorted(2,:)*meter2pixel;
        dx2=x1-x0;
        dy2=y1-y0;
        quiverscale=2;
        quiver(x0,y0,dx2,dy2,quiverscale);
        plot(x1,y1,cfg2);
        title('s-t plane (virtual screen points)');
        
        %% render light fields
        figure(21);
        hold on;
        grid on;
        render_light_field([UV;ST],t_SW,'b-')
        render_light_field([UV_distorted;ST_distorted],t_SW,'r-')
    end
    
    %% append data
    XS_W=[XS_W XS_W0];
    XS_S=[XS_S XS_S0];
    XS_W_distorted=[XS_W_distorted XS_W_distorted0];
    XS_S_distorted=[XS_S_distorted XS_S_distorted0];
    X=[X X0];
    Y=[Y Y0];
    t_E0W=[t_E0W repmat(t_E0W0,1,size(X0,2))];
    X_distorted=[X_distorted X_distorted0];
    G_E=[G_E G_E0];
    G_distorted_E=[G_distorted_E G_distorted_E0];
    
    if 0 % debug
        clf;
        hold on; grid on; axis equal;
        draw_axis(eye(3),zeros(3,1));
        draw_axis(eye(3),t_E0W0);
        plot3v(XS_W0,'r*');
        plot3v(XS_S0,'g*');
        %plotv(XS_W0(1:2,:),'r*');
        %plotv(XS_S0(1:2,:)*meter2pixel,'g*');
        %plotv(XS_S_distorted0(1:2,:)*meter2pixel,'r*');
        keyboard
    end
end
uv_S      = XS_S(1:2,:)          *meter2pixel; % centered pixels in S
uv_dist_S = XS_S_distorted(1:2,:)*meter2pixel; % with distortion

%% Shift and flip screen points
%h0 = 1024;
%h=0;
%w=0;
uv_S(1,:)=uv_S(1,:)+offsetw;
uv_S(2,:)=uv_S(2,:)+offseth;
uv_dist_S(1,:)=uv_dist_S(1,:)+offsetw;
uv_dist_S(2,:)=uv_dist_S(2,:)+offseth;
uv_S(2,:)=h_hmd-uv_S(2,:)-1;
uv_dist_S(2,:)=h_hmd-uv_dist_S(2,:)-1;

%% Plot training data
progn = @(varargin) varargin{end};
myclf = @() progn(...
    evalc('clf'), ...
    evalc('hold on'), ...
    evalc('grid on'), ...
    evalc('axis equal') );

figure(22);
myclf();
plotv(uv_S,'b.');
plotv(uv_dist_S,'r.');
p2d_E0 = uv_S;
duv=uv_dist_S-uv_S;
quiverscale=1;
quiver(p2d_E0(1,:),p2d_E0(2,:),duv(1,:),duv(2,:),quiverscale);
%keyboard;
%return;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xtrain = UVST;
ytrain = UVST_distorted;
%% Regress optical element's distortion
opt='static_xscale';
KR_model_optical=kernel_regression_nDmD(xtrain,ytrain,opt); %(x,y)
KR_model_optical.DISTORT_UVST = DISTORT_UVST;
save 'KR_model_optical' KR_model_optical; %% Save regression parameters
%draw_test_result (xtrain,xtrain, meter2pixel,w,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TODO CHECK the artificial data generation process

%% Create artificial test data
%%% C:\Users\Yuta\Copy\20141214_EyeHMDCalibDataset_ST60_PS\DisplayCalibrationToolbox\data\5\
R_E0W = [
   0.996226308432124   0.031066683714960  -0.081043220263370
  -0.028404458041297   0.999024955127270   0.033798310554464
   0.082014200911078  -0.031368777405414   0.996137375392073
   ];
t_E0W = [
  -0.029793673250077
  -0.078804286612701
   0.118703726732110
   ];
R_WE = R_E0W';
t_WE = -R_E0W'*t_E0W;

grid_step=50;
%[u,v] = meshgrid(1:grid_step:1280,1:grid_step:1024);%%% IEEE VR video
[u,v] = meshgrid(300:grid_step:1080,150:grid_step:800);%%%


%t_EW=-R_WE'*t_EW_dummy;
uv=[u(:) v(:)]';
xyz_gt = R_WE'*(K_E\[uv;ones(1,size(uv,2))]) + repmat(-R_WE'*t_WE,1,size(uv,2));
Y0=repmat(-R_WE'*t_WE,1,size(xyz_gt,2)); % t_EW
%Y0=repmat(t_EW_dummy,1,size(xyz_gt,2));

% dataset for uvst
uv_S2(1,:)=(uv(1,:)-offsetw)/meter2pixel;
uv_S2(2,:)=(uv(2,:)-offseth)/meter2pixel;
X0 = [uv_S2;zeros(1,length(uv))];
X0 = R_SW*X0 + repmat(t_SW,1,size(X0,2));
[XS_W0,XS_S] = intersetRayWithPlane(X0,Y0,t_SW,   R_SW);
[~,UV_test]  = intersetRayWithPlane(X0,Y0,t_SW_z0,R_SW);

xtest = [UV_test(1:2,:);XS_S(1:2,:)];
%%%

yte = kernel_regression_apply(xtest,KR_model_optical);


%% render light fields
figure(34);
clf
hold on; axis equal; grid on; box on;
plot(yte(1,:),yte(2,:),'.');
plot(xtest(1,:),xtest(2,:),'.');
    draw_quivers(xtest,yte);
    legend({'Original','Distorted'});
%
%keyboard

end

%% a function renders s-u plane, t-v plane, and s-t plane
function draw_st_plane (UVST1,UVST2, meter2pixel,w,h,DRAW_CONTOUR,u)
if nargin<6
    DRAW_CONTOUR = false;
end

quiverscale =2;
cfg1 = 'b.';% undistorted
cfg2 = 'r.';% distorted

UV = UVST1(1:2,:);
ST = UVST1(3:4,:);
UV_distorted = UVST2(1:2,:);
ST_distorted = UVST2(3:4,:);
if 0
    ymin = 200;
    %ymax = 1000;
    ymax = 900;
    xmin = 0;
    xmax = 1100;
    xy=ST;
    xy_distorted=ST_distorted;
else
    ymin = 420;
    ymax = 560;
    xmin = 620;
    xmax = 750;
    xy=UV;
    xy_distorted=UV_distorted;
end
%% ST plane
set(gca,'YDir','reverse');
hold on;grid on; axis equal;

xlim([xmin xmax]);ylim([ymin ymax]);

% 1st ST-plane
x0 = xy(1,:)*meter2pixel;
y0 = xy(2,:)*meter2pixel;
x = x0 + w/2-0.5;
y = h - (y0 + h/2-0.5);
plot(x, y, cfg1);
% 2nd ST-plane
x0 = xy_distorted(1,:)*meter2pixel;
y0 = xy_distorted(2,:)*meter2pixel;
x2 = x0 + w/2-0.5;
y2 = h - (y0 + h/2-0.5);
% optical flow between the 1st and 2nd ST planes
px=x2-x;
py=y2-y;
quiverc(x,y,px,py,quiverscale);
plot(x2, y2, cfg2, 'MarkerSize',3);
%title('s-t plane (virtual screen points)');
% draw the contours of the length of the optical flows
if DRAW_CONTOUR
    z=sqrt(sum(px.*px+py.*py,1));
    xx=reshape(x,size(u));yy=reshape(y,size(u));zz=reshape(z,size(u));
    contour(xx,yy,zz);
end
%print(gcf,'filename.pdf','-dpdf','-r0')
return

end
%% a function renders s-u plane, t-v plane, and s-t plane
function draw_test_result (UVST1,UVST2, meter2pixel,w,h,DRAW_CONTOUR,u)
if nargin<6
    DRAW_CONTOUR = false;
end

clf
fig_width =1600;
fig_height=500;
quiverscale =2;
cfg1 = 'b.';
cfg2 = 'r.';
ymin = 200;
%ymax = 1000;
ymax = 900;
xmin = 0;
xmax = 1100;

screen_size=get(gcf,'Position');
set(gcf,'Position',[screen_size(1) screen_size(2) fig_width fig_height]);
subplot_row = 1;
subplot_col = 3;
UV = UVST1(1:2,:);
ST = UVST1(3:4,:);
UV_distorted = UVST2(1:2,:);
ST_distorted = UVST2(3:4,:);

%% U-S plane
subplot(subplot_row,subplot_col,1);
hold on;grid on; axis equal;
plot(UV(1,:),ST(1,:),cfg1);
plot(UV_distorted(1,:),ST_distorted(1,:),cfg2);
legend('Ideal','Distorted');
title('u-s plane')

%% V-T plane
subplot(subplot_row,subplot_col,2);
hold on;grid on; %axis equal;
plot(UV(2,:),ST(2,:),cfg1);
plot(UV_distorted(2,:),ST_distorted(2,:),cfg2);
title('v-t plane');

%% ST plane
subplot(subplot_row,subplot_col,3);
set(gca,'YDir','reverse');
hold on;grid on; axis equal;
xlim([xmin xmax]);ylim([ymin ymax]);

% 1st ST-plane
x0 = ST(1,:)*meter2pixel;
y0 = ST(2,:)*meter2pixel;
x = x0 + w/2-0.5; % OpenCV convention: image origin is (0,0) 
y = h - (y0 + h/2-0.5);
% 2nd ST-plane
x0 = ST_distorted(1,:)*meter2pixel;
y0 = ST_distorted(2,:)*meter2pixel;
x2 = x0 + w/2-0.5;
y2 = h - (y0 + h/2-0.5);
% optical flow between the 1st and 2nd ST planes
px=x2-x;
py=y2-y;

if 1 % IEEE VR, pre-distored grid image for ST60
    w=1280;
    h=1024;
    image_org = uint8(zeros(h,w,3));
    image_est = uint8(zeros(h,w,3));
    xidx=int32(x)+1;
    yidx=int32(y)+1;
    x2idx=int32(x2)+1;% distorted
    y2idx=int32(y2)+1;
    idx= (w>=x2idx) & (x2idx>0) & (h>=y2idx) & (y2idx>0);
    x2idx=x2idx(idx);
    y2idx=y2idx(idx);
    xidx=xidx(idx);
    yidx=yidx(idx);
    %for k=1:length(x2)
        image_org(yidx, xidx, 2)=255;
        image_org(y2idx,x2idx,1)=255;
        %image_est(y2idx,x2idx,1)=255;
    %end
    imshow(image_org)
    %imshow(image_est)
    %%keyboard
else
    plot(x, y, cfg1);
    plot(x2, y2, cfg2);
    quiverc(x,y,px,py,quiverscale);
    title('s-t plane (virtual screen points)');
    % draw the contours of the length of the optical flows
    if DRAW_CONTOUR
        z=sqrt(sum(px.*px+py.*py,1));
        xx=reshape(x,size(u));yy=reshape(y,size(u));zz=reshape(z,size(u));
        contour(xx,yy,zz);
    end
    %print(gcf,'filename.pdf','-dpdf','-r0')
end
return
end


%%
function draw_axis(R,t,s,color,kLineWidth)
if nargin<=0
    R=eye(3);
    t=[0 0 0]';
    s=0.5;
    color=eye(3);
    kLineWidth=1;
elseif nargin<=1
    t=[0 0 0]';
    s=0.5;
    color=eye(3);
    kLineWidth=1;
elseif nargin<=2
    s=0.5;
    color=eye(3);
    kLineWidth=1;
elseif nargin<=3
    color=eye(3);
    kLineWidth=1;
elseif nargin<=4
    kLineWidth=1;
    
end
o=[0 0 0]'+t;
x=R*[s 0 0]'+t;
y=R*[0 s 0]'+t;
z=R*[0 0 s]'+t;
tmp=[o x]';
cfg='-';
plot3(tmp(:,1),tmp(:,2),tmp(:,3),cfg,'LineWidth',kLineWidth,'color',color(1,:));
tmp=[o y]';
plot3(tmp(:,1),tmp(:,2),tmp(:,3),cfg,'LineWidth',kLineWidth,'color',color(2,:));
tmp=[o z]';
plot3(tmp(:,1),tmp(:,2),tmp(:,3),cfg,'LineWidth',kLineWidth,'color',color(3,:));
end

%% Project and plot 3D points xyz by applying the projection matrix P
function xy=projection(P,xyz)
xyzw=[xyz; ones(1,length(xyz))];
xyz = P*xyzw;

xy = [xyz(1,:)./xyz(3,:); xyz(2,:)./xyz(3,:)];
end

%%
function x_n = normalize(x,m,std)
assert( size(x,1)==size(m,1) && size(m,1)==size(std,1), 'input should be x: DxN, m:Dx1, and s:Dx1' )
N=size(x,2);
x_n = (x - repmat(m,1,N))./repmat(std,1,N);
end

%%
function x = unnormalize(x_n,m,ssqrt)
assert( size(x_n,1)==size(m,1) && size(m,1)==size(ssqrt,1), 'input should be x: DxN, m:Dx1, and s:Dx1' )
N=size(x_n,2);
x = ( x_n.*repmat(ssqrt,1,N) ) + repmat(m,1,N);
end

%%
function R = skew(w)
%SKEW  generates a skew-symmetric matrix given a vector w
%
%	R = SKEW(w)
%
% See also: ROTAXIS, SKEWEXP, SKEWCOORDS.

% $Id: skew.m,v 1.1 2009-03-17 16:40:18 bradleyk Exp $
% Copyright (C) 2005, by Brad Kratochvil

if 3 ~= size(w,1),
    error('SCREWS:skew','vector must be 3x1')
end

if isnumeric(w),
    R = zeros(3,3);
end

R(1,2) = -w(3);
R(1,3) =  w(2);
R(2,3) = -w(1);

R(2,1) =  w(3);
R(3,1) = -w(2);
R(3,2) =  w(1);

%   R(1,1) = 0;
%   R(2,2) = 0;
%   R(3,3) = 0;

end

%% Render light field
function render_light_field(UVST,t_SW,cfg)

        UV = UVST(1:2,:);
        ST = UVST(3:4,:);
        x1=[ST(1,:);ST(2,:);t_SW(3)*ones(1,length(UV))];
        x2=[UV(1,:);UV(2,:);zeros(1,length(UV))];
        N=length(UV);
        for j=1:1:N
            tmp=[x1(:,j) x2(:,j)]';
            plot3(tmp(:,1),tmp(:,2),tmp(:,3),cfg);
        end
end

%% Rodrigues formula
function R = rodrig(k,fi)
% This is just to make it easyer to read!
x = k(1);
y = k(2);
z = k(3);

% Create a 3x3 zero matrix
R = zeros(3,3);
% We use the formual for rotationg matrix about a unit vector k

R(1,1) = cos(fi)+x^2*(1-cos(fi));
R(1,2) = x*y*(1-cos(fi))-z*sin(fi);
R(1,3) = x*z*(1-cos(fi))+y*sin(fi);

R(2,1) = y*x*(1-cos(fi))+z*sin(fi);
R(2,2) = cos(fi)+y^2*(1-cos(fi));
R(2,3) = y*z*(1-cos(fi))-x*sin(fi);

R(3,1) = z*x*(1-cos(fi))-y*sin(fi);
R(3,2) = z*y*(1-cos(fi))+x*sin(fi);
R(3,3) = cos(fi)+z^2*(1-cos(fi));
end