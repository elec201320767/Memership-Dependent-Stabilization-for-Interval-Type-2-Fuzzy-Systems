%% Sparse multiplication matrix


function sparse = E4x4(e13x3,e23x3,e33x3,e43x3,i)
    if i == 1
        sparse = e13x3;
    elseif i == 2
            sparse = e23x3;
    elseif i == 3
        sparse = e33x3;
    else
        sparse = e43x3;
    end
end