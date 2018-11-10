function [u, f] = tv_deblur_noncyclic(u, k, lambda, mu, tmax)
% ---- Description:
% Use total variation model to deblur and denoise an image. See
% doc\report.pdf for further information. Here we assume the convolutional
% operator is cyclic
% ---- Inputs:
%   u: the origin image
%   k: the blur kernel
%   lambda: the coefficent for the regularization term
%   mu: the coefficent for the second order penalizing ||d - Wu||_2^2,
%   tmax: the maximum allowable iteration
% ---- Outputs:
%   u: The recovered image
%   f: The blurred image

f = imfilter(u, k);
sigma1 = max(max(max(abs(u))))/100;
f = f + sigma1*randn(size(f));
u = f;

% gradient operator
grad_u1 = imfilter(f, [1,-1]);
grad_u2 = imfilter(f, [1;-1]);

% initialization
p1 = zeros(size(grad_u1,1), size(grad_u1,2));
p2 = zeros(size(grad_u2,1), size(grad_u2,2));
tau = lambda/mu;
ktf = imfilter(f(end:-1:1, end:-1:1, :),k);
ktf = ktf(end:-1:1, end:-1:1,:); % ktf is a variable represents A^t*b in the algorithm
op= @(x) operate(x,mu, k, [1,-1], [1;-1]); % corresponding to the operator (A^tA + mu(W1^t W1 + W2^t W2))

% solve the problem using admm
for t = 1:tmax
    v1 = grad_u1 + 1/mu*p1;
    d1 = sign(v1).*max(abs(v1) - tau,0);
    v2 = grad_u2 + 1/mu*p2;
    d2 = sign(v2).*max(abs(v2) - tau,0);
    temp1 = d1 - 1/mu*p1;
    temp2 = d2 - 1/mu*p2;
    Wtemp1 = imfilter(temp1(end:-1:1, end:-1:1, :),[1,-1]);
    Wtemp1 = Wtemp1(end:-1:1, end:-1:1,:); % W1^t*temp1
    Wtemp2 = imfilter(temp2(end:-1:1, end:-1:1, :),[1;-1]);
    Wtemp2 = Wtemp2(end:-1:1, end:-1:1,:);
    b = ktf + mu*Wtemp1 + mu*Wtemp2;  % W2^t*temp2
    u = cg(op, b, u);
    grad_u1 = imfilter(u, [1,-1]);
    grad_u2 = imfilter(u, [1;-1]);
    p1 = p1 + mu*(grad_u1 - d1);
    p2 = p2 + mu*(grad_u2 - d2);
    %disp(sum(sum(sum(abs(grad_u1 - d1)))) + sum(sum(sum(abs(grad_u2 - d2)))))
end
end

