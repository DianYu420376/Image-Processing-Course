function [u, f] = tv_deblur_cyclic(u, k, lambda, mu, tmax)
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


fu = fft2(u);
fk = psf2otf(k, [size(u,1),size(u,2)]);
ff = fu.*fk;
f = real(ifft2(ff));
sigma1 = max(max(max(abs(u))))/100;
f = f + sigma1*randn(size(f));
ff = fft2(f);

% gradient operator
grad1 = [1, -1];
grad2 = [1; -1];
fgrad1 = psf2otf(grad1,[size(u,1),size(u,2)]);
fgrad2 = psf2otf(grad2,[size(u,1),size(u,2)]);
grad1_u = real(ifft2(fgrad1.*ff));
grad2_u = real(ifft2(fgrad2.*ff));

% initialization
p1 = zeros(size(u,1), size(u,2));
p2 = zeros(size(u,1), size(u,2));
tau = lambda/mu;
precomputed = conj(fk).*fk+ mu*(conj(fgrad1).*fgrad1 + conj(fgrad2).*fgrad2);

% solve the problem using admm
for t = 1:tmax
    v1 = grad1_u + 1/mu*p1;
    d1 = sign(v1).*max(abs(v1) - tau,0);
    v2 = grad2_u + 1/mu*p2;
    d2 = sign(v2).*max(abs(v2) - tau,0);
    ftemp1 = fft2(d1 - 1/mu*p1);
    ftemp2 = fft2(d2 - 1/mu*p2);
    fu = (conj(fk).*ff + mu* (conj(fgrad1).* ftemp1 + conj(fgrad2).*ftemp2))./precomputed;
    grad1_u = real(ifft2(fgrad1.*fu));
    grad2_u = real(ifft2(fgrad2.*fu));
    p1 = p1 + mu*(grad1_u - d1);
    p2 = p2 + mu*(grad2_u - d2);
    
    % calculating the primal and dual value and view the duality gap
    primal = 1/2*sum(sum(sum(abs(fk.*fu - ff).^2))) + lambda*sum(sum(sum(abs(d1)))) + ...
        lambda*sum(sum(sum(abs(d2))));
    fp1 = (1./conj(fk)).*(conj(fgrad1)).*fft2(p1);
    fp2 = (1./conj(fk)).*(conj(fgrad2)).*fft2(p2);
    fp = fp1 + fp2;
    dual =   -1/2*sum(sum(sum(abs(fp).^2))) +1/2*sum(sum(sum(conj(ff).*(fp) + conj(fp).*ff)));
    %disp(primal - dual)
    %disp(sum(sum(sum(abs(grad1_u - d1)))) + sum(sum(sum(abs(grad2_u - d2)))))
end
u = real(ifft2(fu));
end

