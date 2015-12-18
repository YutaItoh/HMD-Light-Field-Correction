
function h=plotv(varargin)
assert( nargin >= 1, 'Usage: plotv(X,[option])');
v = varargin{1};
assert(size(v,1)==2, 'X has to be a 2xN matrix');
if nargin >= 2
    h=plot(v(1,:),v(2,:),varargin{2:end});
else
    h=plot(v(1,:),v(2,:));
end
end