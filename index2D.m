function index = index2D(M)

%    M = 2
%
%    #   index                Norm
%      +------------------------
%    1 |  0      0      0        0
%      |
%    2 |  0      0      1        1
%    3 |  0      1      0        1
%    4 |  1      0      0        1
%      |
%    5 |  0      0      2        2
%    6 |  0      1      1        2
%    7 |  1      0      1        2
%    8 |  0      2      0        2
%    9 |  1      1      0        2
%   10 |  2      0      0        2
%      


% index = zeros(N,3);

i = 1;
idxnorm = M;

index(i,:) = [0,0,idxnorm];
t3 = index(i,3);
t2 = index(i,2);
t1 = index(i,1);
i = i+1;

while t1<=idxnorm

    if (t3+t1 == idxnorm)
        t3 = t3-1;
        t1 = 0;
        t2 = idxnorm - t3;
    else
        t2 = t2-1;
        t1 = idxnorm - t3 - t2;
    end

    if t3 == -1
        break
    end

    index(i,:) = [t1,t2,t3];
    i = i+1;
end