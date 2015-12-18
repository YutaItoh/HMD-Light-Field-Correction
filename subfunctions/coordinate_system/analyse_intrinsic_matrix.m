function [convention_enum ,K2]= analyse_intrinsic_matrix(K)
    display( strcat('The camera coordinate system is assumed to be RIGHT-HANDED.') );
    assert(size(K,1)==3 && size(K,2)==3, 'Input must be 3x3 matrix');
    assert( K(1,1)~=0, 'K(1,1) should be non-zero');
    assert( K(2,2)~=0, 'K(2,2) should be non-zero');
    assert( K(3,3)~=0, 'K(3,3) should be non-zero');

    R = eye(3);
    is_forwarding_z = false;
    is_topleft = true;
    
    %% Warnings
    if K(1,1) < 0
        K = -K;
        warning('K(1,1) is negative, -K is used to make sure K(1,1)>0');
    end
    if abs(K(3,3))~=1 
        msg = strcat('abs(K(3,3)) = ', {' '}, num2str(K(3,3), 17) ,' is not 1, is this intended?');
        warning(char(msg));
    end
    if sign(K(1,3))*sign(K(3,3)) <0
        display( K(1,3)*K(3,3) );
        msg = 'The x value of the image center will be negative, is this intented?';
        warning(msg);
    end
    if sign(K(2,3))*sign(K(3,3)) <0
        display( K(2,3)*K(3,3) );
        msg = 'The y value of the image center will be negative, is this intented?';
        warning(msg);
    end
    % Z-axis direction
    if K(3,3)>0
        is_forwarding_z = true;
         R = [1 0 0;0 -1 0;0 0 -1];
    end    
    % Image plane origin
    if ( K(3,3)<0  && K(2,2)>0 ) || ( K(3,3)>0  && K(2,2)<0 )
        is_topleft = false;
    end
    
    display(K);
    
    if is_forwarding_z==true  && is_topleft==true
        display( 'K follows OpenCV convention.' );
        convention_enum = CoordinateConventions.OpenCV;
    end
    if is_forwarding_z==false && is_topleft==false
        display( 'K follows Ubitrack convention.' );
        convention_enum = CoordinateConventions.Ubitrack;
    end
    if is_forwarding_z==false && is_topleft==true
        display( 'K follows Our (OpenGL) convention.' );
        convention_enum = CoordinateConventions.OpenGL;
    end
    if is_forwarding_z==true && is_topleft==false
        display( 'K follows uncommon convention.' );
        convention_enum = CoordinateConventions.Uncommon;
    end
    
    %% Visualize the camera
%     clf;
%     hold on; grid on;axis equal;
%     
%     len=1;
%     draw_axis(eye(3),ones(3,1),len,eye(3),5); % Reference coordinate system
%     draw_axis(R,zeros(3,1)); % Camera coordinate system
%     w=6.4;
%     h=4.8;
%     a=1;
%     draw_screen(w,h,a,eye(3),[-w/2;-h/2;-1],'b',is_topleft ); % Image plane

    K2=K;
end

    