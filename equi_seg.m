function X = equi_seg(N)


X=[];
for i =0:N-1
    for j=0:N-i-1
        X(end+1,:) = [i/N,j/N];
        X(end+1,:) = [(i+1)/N,j/N];
        X(end+1,:) = [i/N,j/N];
        X(end+1,:) = [i/N,(j+1)/N];
    end
end

for k = 0:N
        for i = 0:k-1
            X(end+1,:) = [i/N,(k-i)/N];
            X(end+1,:) = [(i+1)/N,(k-i-1)/N];
        end
end