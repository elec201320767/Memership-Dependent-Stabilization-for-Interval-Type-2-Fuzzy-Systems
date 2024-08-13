
%%
%% Vertices3D calculation for T1 FS
function vertices = ExtremaCombination(p,Extrema)
%% Combination 2^(p-1) of extrema
% p > 2
comb = 1: 2*(p-1);
comb = nchoosek(comb,p-1);
temp = mod(comb,p-1);
Sizecomb = size(temp);
index = [];
for i=Sizecomb(1):-1:1
    for j=1:p-1
        for k=1:p-1
            if k == j
                continue
            end
            if temp(i,j) == temp(i,k)
                index = [index, i];
                continue
            end
        end
    end
end
index = unique(index);
Sizeindex = size(index);
for i=Sizecomb(1):-1:1
    for j = 1:Sizeindex(2)
        if i==index(j)
            comb(i,:) = [];
        end
    end
end
%rearange comb
temp = mod(comb,p-1);
Sizecomb = size(temp);
for i = 1:Sizecomb(1)
    for j = 1:Sizecomb(2)
        if temp(i,j) == 0
            temp(i,j) = p-1;
        end
    end
end
for i = 1:Sizecomb(1)
    for j = 1:Sizecomb(2)
        Comb(i,temp(i,j)) = comb(i,j);
    end
end
% comb is [2^p x p] size

for i = 1:2^(p-1)
    for j = 1:p-1
       ExtremeRect(i,j) = Extrema(1,Comb(i,j)) ;
    end
end
vertices  = ExtremeRect;
end
