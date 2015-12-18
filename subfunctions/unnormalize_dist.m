function x = unnormalize_dist(x_n,m,ssqrt)
assert( size(x_n,1)==size(m,1) && size(m,1)==size(ssqrt,1), 'input should be x: DxN, m:Dx1, and s:Dx1' )
N=size(x_n,2);
x = ( x_n.*repmat(ssqrt,1,N) ) + repmat(m,1,N);
end
