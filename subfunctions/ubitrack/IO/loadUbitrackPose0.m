% Convert "Print Sink (Quaternion)" output format to a rotation matrix
function [R, t, q]  = loadUbitrackPose0(file)
    function R = ubitrackQuat2Mat(q)
        if(size(q,1) == 1)
            q=q';
        end
        I=[0 0 -1 0; 0 -1 0 0 ; 1 0 0 0; 0 0 0 1];
        q;
        R = diag([-1 -1 1])*quat2dcm( (I*q)')';
        %R = quat2dcm( q');
    end
tmp=openUbitrack6DPoseCalib(file);
tstamp =tmp(1);
q = tmp(2:5);
R = ubitrackQuat2Mat(q);
t = tmp(end-2:end);
%t=-R'*t;
%R=R' %%% Needed????
end