
%% 2022.04  kyungsoo@postech.ac.kr
%% Vertices3D calculation for T2 Fuzzy System
function verticeSet = Type2FuzzyVerticesBound(MF,mfsum)
%% define  a membership function(MF)
%  MF has a form in [Lower MF ; Upper MF] 
%% Vertices of Extreme box
SizeMF = size(MF);
p = SizeMF(1,1)/2;      % # of rules
comb = 1: 2*p;          
 % to obtain vertices of exrtreme rectangular box
  
for i = 1:p
    a(1,i) = min(MF(i,:));
    b(1,i) = max(MF(p+i,:));
end

comb = nchoosek(comb,p);
temp = mod(comb,p);
Sizecomb = size(temp);
index = [];
for i=Sizecomb(1):-1:1
    for j=1:p
        for k=1:p
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
temp = mod(comb,p);
Sizecomb = size(temp);
for i = 1:Sizecomb(1)
    for j = 1:Sizecomb(2)
        if temp(i,j) == 0
            temp(i,j) = p;
        end
    end
end
for i = 1:Sizecomb(1)
    for j = 1:Sizecomb(2)
        Comb(i,temp(i,j)) = comb(i,j);
    end
end

% comb is [2^p x p] size
Extrema = [a b];

for i = 1:2^p
    for j = 1:p
       ExtremeRect(i,j) = Extrema(1,Comb(i,j)) ;
    end
end   
% ExtremeRect has extrema vertices with [2^p x p] size
%% Vertices Algorithm
v = 0;
vertices = zeros(1,p);
  for i= 1:p
    ExtremaComb  = Extrema;
    ExtremaComb(i+p) = [];
    ExtremaComb(i) = [];
    RectVertices = ExtremaCombination(p,ExtremaComb);     % 2^(p-1) combinations of extrema of MF for j \neq i, j = 1,...,p
    for j = 1:2^(p-1)
        convSpace = mfsum - sum( RectVertices(j,:));
        if a(i) <= convSpace & convSpace <= b(i) % vertices condition
            v = v + 1;
            temp = 0;
            for k = 1:p
                if k == i
                    vertices(v,k) = convSpace;
                    temp = 1;
                else
                    vertices(v,k) = RectVertices(j,k-temp);
                end
            end
        end
    end
  end 
 verticeSet = round(vertices,4);
 verticeSet = unique(verticeSet,'rows');
 verticeSet = verticeSet';
 
 ExtremeRect = ExtremeBox(MF);
 icept = sum(ExtremeRect');
 
 licept = find(1>icept);
 lExtremeRect = ExtremeRect(licept,:);
 uicept = find(1<icept);
 uExtremeRect = ExtremeRect(uicept,:); 
 if mfsum<1 % LMF
   lindex = find(mfsum < sum(lExtremeRect'));
   verticeSet = [verticeSet lExtremeRect(lindex,:)'];
 else       % UMF
   uindex = find(mfsum > sum(uExtremeRect'));
   verticeSet = [verticeSet uExtremeRect(uindex,:)']; 
 end
 

  
 
 
 
end
