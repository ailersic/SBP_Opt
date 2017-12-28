function [c, ceq] = nonlin_con(Hsol, d, d_)
    c = [-subs(diag(Hsol), d, d_);
         -d_
          d_ - ones(length(d_), 1)];
    for i=1:length(d_)-1
        c = [c; d_(i) - d_(i+1)];
    end
    c = double(c);
    ceq = [];
end