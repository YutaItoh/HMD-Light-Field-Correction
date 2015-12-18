function [X_W,Y,X_distorted, t_E0W, t_WE0, G_E, G_distorted_E]=loadDataSet2(dir,K_E, misc)

USE_NEW_DATAFORMAT=misc.USE_NEW_DATAFORMAT;
if USE_NEW_DATAFORMAT
    dir
    filename = strcat(dir,'3dPointsOnGrid.txt');
    p3D=openUbitrack3DPositionList(filename);
    filename = strcat(dir,'2dPointsOnViewPoint.txt');
    p2D_E = openUbitrack2DPositionList(filename); % marker 2D points in undistorted image
    p2D_E(2,:) = misc.h-1-p2D_E(2,:);
    filename = strcat(dir,'2dPointsOnViewPointThroughHMD.txt');
    p2D_distorted_E=openUbitrack2DPositionList(filename); % marker 2D points through HMD in undistorted image
    p2D_distorted_E(2,:) = misc.h-1-p2D_distorted_E(2,:);
    filename = strcat(dir,'P_HMDCam2IDS.txt');
    [R_E0W, t_E0W]=loadUbitrackPose0(filename)
    t_WE0 = -R_E0W'*t_E0W;
    filename = strcat(dir,'P_IDS2Marker.txt');
    [R_ME0, t_ME0]=loadUbitrackPose0(filename);
else
    filename = strcat(dir,'../p3d.txt');
    p3D=openUbitrack3DPositionList(filename);
    filename = strcat(dir,'p2d_E0.txt');
    p2D_E = openUbitrack2DPositionList(filename); % marker 2D points in undistorted image
    p2D_E(2,:) = misc.h-1-p2D_E(2,:);
    filename = strcat(dir,'p2d_E1.txt');
    p2D_distorted_E=openUbitrack2DPositionList(filename); % marker 2D points through HMD in undistorted image
    p2D_distorted_E(2,:) = misc.h-1-p2D_distorted_E(2,:);
    filename = strcat(dir,'pose_E0W_error.txt');
    [R_E0W, t_E0W]=loadUbitrackPose0(filename);
    t_WE0 = -R_E0W'*t_E0W;
    filename = strcat(dir,'pose_ME0_error.txt');
    [R_ME0, t_ME0]=loadUbitrackPose0(filename);
end

P_E0W=[ [R_E0W t_E0W]; 0 0 0 1];
P_ME0=[ [R_ME0 t_ME0]; 0 0 0 1];
P_MW = P_E0W*P_ME0;
R_MW = P_MW(1:3,1:3);
t_MW = P_MW(1:3,4);

N=size(p3D,2);
X_W = R_MW*p3D+repmat(t_MW,1,N); % grid points
Y=repmat(t_E0W,1,N);
G_distorted_E = K_E\[p2D_distorted_E;ones(1,N)];
X_distorted=R_E0W*G_distorted_E + repmat(t_E0W,1,N);
%X_distorted=X_distorted/max(X_distorted(:))*0.3;


        
G_E = R_ME0*p3D+repmat(t_ME0,1,N); % grid points
uv_E=project(K_E,G_E);
% uv_E and p2D_E should be almost idential


if 0
    figure(misc.fig_id);
    clf
    hold on; grid on; axis equal;
    xlim([misc.xlim_min,misc.xlim_max]);
    ylim([misc.ylim_min,misc.ylim_max]);
    plot(p2D_E(1,:),p2D_E(2,:),'ro')
    plot(p2D_distorted_E(1,:),p2D_distorted_E(2,:),'b.')
    duv=p2D_distorted_E-p2D_E;
    quiverscale=1;
    quiver(p2D_E(1,:),p2D_E(2,:),duv(1,:),duv(2,:),quiverscale);
    legend('Original','Distorted');
    title('Training data')
    keyboard
else
    if USE_NEW_DATAFORMAT
        img_filename = strcat(dir,'UserViewUndistorted.png');
    else
        img_filename = strcat(dir,'image_eye_undist_before.png');
    end
    %img_before = imread(img_filename);
    
    
    if 0 %debug
    figure(21); clf;
    imshow(img_before), hold on;
    kLineWidth  = 2;
    kMarkerSize = 2;
    plot([0 100],[0 100],'go-','LineWidth', kLineWidth,'MarkerSize',kMarkerSize)
    plot(p2D_E(1,:),p2D_E(2,:),'go','LineWidth', kLineWidth,'MarkerSize',kMarkerSize)
    plot(uv_E(1,:),uv_E(2,:),'bo','LineWidth', kLineWidth,'MarkerSize',kMarkerSize)
    plot(p2D_distorted_E(1,:),p2D_distorted_E(2,:),'ro','LineWidth', kLineWidth,'MarkerSize',kMarkerSize)
    end
    
    if 0 %debug
        figure(22);
        dif=uv_E - p2D_E;
        clf
        subplot(1,2,1);
        hold on; grid on; axis equal;
        plot(dif(1,:),dif(2,:),'*r')
        dif=p2D_distorted_E - p2D_E;
        subplot(1,2,2);
        hold on; grid on; axis equal;
        plot(dif(1,:),dif(2,:),'*r')
    end
%keyboard
end
return

%dir = '.\';
    if USE_NEW_DATAFORMAT
        display_param_linear = load('C:/Users/Yuta/Dropbox/work/Projects/20140530_EyeHMDCalibDataset/calib_hmd_optics/display_param_linear.mat');

        %display_param_linear = load(strcat(dir,'../display_param'));
        display_param_nonlin = display_param_linear;
    else
display_param_nonlin = load(strcat(dir,'../display_param_nonlinear'));
display_param_linear = load(strcat(dir,'../display_param_linear'));
    end

if 1 % convert coordinate system
    display_param_nonlin.R_WS(2:3,:) = -display_param_nonlin.R_WS(2:3,:);
    display_param_nonlin.t_WS(2:3)   = -display_param_nonlin.t_WS(2:3);
    display_param_linear.R_WS(2:3,:) = -display_param_linear.R_WS(2:3,:);
    display_param_linear.t_WS(2:3)   = -display_param_linear.t_WS(2:3);
end

R_WS = display_param_linear.R_WS;
t_WS = display_param_linear.t_WS;
meter2pixel = display_param_linear.alpha;
%keyboard

end