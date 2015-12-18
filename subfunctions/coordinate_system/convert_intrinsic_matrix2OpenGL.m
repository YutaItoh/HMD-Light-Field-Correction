%% Convert intrinsic so that it compatible to our (OpenGL) convention
    function K_new = convert_intrinsic_matrix2OpenGL(K,h)
        [source_convention,K_new] = analyse_intrinsic_matrix(K);
%     if is_forwarding_z==true  && is_topleft==true
%         display( 'K follows OpenCV convention.' );
%         convention_enum = CoordinateConventions.OpenCV;
%     end
%     if is_forwarding_z==false && is_topleft==false
%         display( 'K follows Ubitrack convention.' );
%         convention_enum = CoordinateConventions.Ubitrack;
%     end
%     if is_forwarding_z==false && is_topleft==true
%         display( 'K follows Our (OpenGL) convention.' );
%         convention_enum = CoordinateConventions.OpenGL;
%     end
%     if is_forwarding_z==true && is_topleft==false
%         display( 'K follows uncommon convention.' );
%         convention_enum = CoordinateConventions.Uncommon;
%     end
    if     source_convention == CoordinateConventions.OpenCV;
        K_new(:,2:3) = -K_new(:,2:3);
    elseif source_convention == CoordinateConventions.Ubitrack;
        % flip the image center
        K_new(2,3) = sign(K_new(3,3))*(h-abs(K_new(2,3)) -1 );
        K_new(2,2)   = -K_new(2,2);
    elseif source_convention == CoordinateConventions.OpenGL;
        % do nothing
    elseif source_convention == CoordinateConventions.Uncommon;
        K_new(2,3) = sign(K_new(3,3))*(h-sign(K_new(3,3))*K_new(2,3) -1 );
        K_new(:,3)   = -K_new(:,3);
    end

    analyse_intrinsic_matrix(K_new);
end

    