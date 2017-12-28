function d0 = gen_d0(m, i, nd)
    d0 = [];
    for j=i:m-nd+1
        i1 = j/(m+1);
        if nd > 1
            in = gen_d0(m, j+1, nd-1);
        else
            in = [];
        end
        d0 = [d0, [i1*ones(1, max(size(in, 2), 1)); in]];
    end
end
