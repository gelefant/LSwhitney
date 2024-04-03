function index = index3D(M)

%    M = 2
%
%    #   index                Norm
%      +------------------------
%    1 |  0      0      0      0        0
%      |
%    2 |  0      0      0      1        1
%    3 |  0      0      1      0        1
%    4 |  0      1      0      0        1
%    5 |  1      0      0      0        1
%      |
%    6 |  0      0      0      2        2
%    7 |  0      0      1      1        2
%    8 |  0      1      0      1        2
%    9 |  1      0      0      1        2
%   10 |  0      0      2      0        2
%   11 |  0      1      1      0        2
%   12 |  1      0      1      0        2
%   13 |  0      2      0      0        2
%   14 |  1      1      0      0        2
%   15 |  2      0      0      0        2
%      


% index = zeros(N,3);

i = 1;
idxnorm = M;

index(i,:) = [0,0,0,idxnorm];
t4 = index(i,4);
t3 = index(i,3);
t2 = index(i,2);
t1 = index(i,1);
i = i+1;

while t1<=idxnorm

    if (t4+t1 == idxnorm)
        t4 = t4-1;
        t1 = 0;
        t2 = 0;
        t3 = idxnorm - t4;
    else
        if (t4+t3+t1 == idxnorm)
            t3 = t3-1;
            t2 = idxnorm - t4 - t3;
            t1 = 0;
        else
        t2 = t2-1;
        t1 = idxnorm - t4 - t3 - t2;
        end
    end

    if t4 == -1
        break
    end

    index(i,:) = [t1,t2,t3,t4];
    i = i+1;
end