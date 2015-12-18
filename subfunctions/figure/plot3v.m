
function h=plot3v(varargin)
assert( nargin >= 1, 'Usage: plot3v(X,[option])');
v = varargin{1};
assert(size(v,1)==3, 'X has to be a 3xN matrix');
if nargin >= 2
    h = plot3(v(1,:),v(2,:),v(3,:),varargin{2:end});
else
    h = plot3(v(1,:),v(2,:),v(3,:));
end

end