function [c, ceq] = nonlin_con(H, d, d_)
    c = [-subs(diag(H), d, d_);
         -d_
          d_ - ones(length(d_), 1)];
    for i=1:length(d_)-1
        c = [c; d_(i) - d_(i+1)];
    end
    c = double(c);
    ceq = [];
end