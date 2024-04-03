function X = equi_seg_int(N)


X=[];
for i =0:N-1
    for j=0:N-i-1
        X(end+1,:) = [i/N,j/N];
        X(end+1,:) = [(i+1)/N,j/N];
        X(end+1,:) = [i/N,j/N];
        X(end+1,:) = [i/N,(j+1)/N];
    end
end

for i = 0:N-1
    X(end+1,:) = [i/N,(N-i)/N];
    X(end+1,:) = [(i+1)/N,(N-i-1)/N];
end


end