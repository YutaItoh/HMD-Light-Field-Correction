function uv = project(K,xyz)
    assert(size(xyz,1) == 3);
    assert(size(K,1) == 3);
    assert(size(K,2) == 3 || size(K,2) == 4);
    if size(K,2) == 4
        xyz=[xyz;ones(1,size(xyz,2))];
    end
    uv=K*xyz;
    uv=[uv(1,:)./uv(3,:); uv(2,:)./uv(3,:)];
end