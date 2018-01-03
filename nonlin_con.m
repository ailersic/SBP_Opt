function [c, ceq] = nonlin_con(H, d, d_, ceqf, desvar)
    c = [-subs(diag(H), d, d_);
         -d_
          d_ - ones(length(d_), 1)];
    for i=1:length(d_)-1
        c = [c; d_(i) - d_(i+1)];
    end
    c = double(c);
    if isempty(ceqf)
        ceq = [];
    else
        ceq = ceqf(d_) - desvar;
    end
end