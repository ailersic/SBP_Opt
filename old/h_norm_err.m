function err = h_norm_err(n, d)
    % x is on domain [0,1]
    if mod(n, 2) == 1
        x = [0; d; 0.5; 1.-flipud(d); 1];
    else
        x = [0; d; 1.-flipud(d); 1];
    end

    [H, Q] = acc_eqn_d(n, d);
    j = n-1;
    e = Q*x.^j - j*H*x.^(j-1);
    
    err = double(e.'*H*e);
end

