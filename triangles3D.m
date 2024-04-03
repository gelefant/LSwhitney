function [Xmin,Xint] = triangles3D(r)

X = [];
Y = [];
for i = 0:r-1
    T = triangles2D(r-i)*(r-i)/r;
    X = [X; [T,repmat(i/r,size(T,1),1)]];
    X = [X; [repmat(i/r,size(T,1),1),T]];
    X = [X; [T(:,1),repmat(i/r,size(T,1),1),T(:,2)]];
    Y = [Y; [T,repmat(i/r,size(T,1),1)]];
end

B = [-1, -1,-1;1,0,0;0,1,0]; b = [1,0,0];

Y = Y*B'+repmat(b,size(Y,1),1);

Xmin = [X;Y];

TT = triangles2D(r);

Xint = [X;[TT,zeros(size(TT,1),1)]*B'+repmat(b,size(TT,1),1)];

end