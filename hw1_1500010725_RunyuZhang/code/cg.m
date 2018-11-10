function x = cg(op, b, x0)
% using conjugate gradient to solve linear equation system
x = x0;
Ax = op(x);
r = b - Ax;
a1 =dot(vec(r), vec(r));
p = r;
Ap = op(r);
a2 = dot(vec(p),vec(Ap));
alpha = a1/a2;
e = norm(vec(b - Ax));
while e > 1e-3
    x = x + alpha*p;
    Ax = op(x);
    r = b - Ax;
    a2 = a1;
    a1 =dot(vec(r), vec(r));
    beta = a1/a2;
    p = r + beta*p;
    Ap = op(p);
    a2 = dot(vec(p),vec(Ap));
    alpha = a1/a2;
    e = norm(vec(b - Ax));
    %disp(e)
end
end

function y = vec(x)
y = reshape(x, [numel(x),1]);
end