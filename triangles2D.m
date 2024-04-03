function [X] = triangles2D(r)

T = [0,0;1,0;0,1];

p = 0:r;
[xp,yp] = meshgrid(p);
xp = xp(:); yp = yp(:);

id = (xp+yp)<=r-1;

xp = xp(id); yp = yp(id);

PPx = repmat(xp'/r,3,1);
PPy = repmat(yp'/r,3,1);
PP = [PPx(:),PPy(:)];
TT = repmat(T/r,size(xp,1),1);

X = TT+PP;

end