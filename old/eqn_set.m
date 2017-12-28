function err = eqn_set(eqns, g)
    g = num2cell(g);
    err = sum(eqns(g{:}).^2);
end