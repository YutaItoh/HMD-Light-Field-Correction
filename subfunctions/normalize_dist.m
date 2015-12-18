
function x_n = normalize_dist(x,m,std)
assert( size(x,1)==size(m,1) && size(m,1)==size(std,1), 'input should be x: DxN, m:Dx1, and s:Dx1' )
N=size(x,2);
x_n = (x - repmat(m,1,N))./repmat(std,1,N);
end