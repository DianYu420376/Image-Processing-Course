function y = operate(x, mu, k, h, v)
% Define the operation 'A^tA + mu*(W1^t*W1 + W2^t*W2)'
% here A corresponds to the convolutional operator k
% W1 corresponds to h
% W2 corresponds to v
y1 = imfilter(x, k);
y1 = imfilter(y1(end:-1:1, end:-1:1, :),k);
y1 = y1(end:-1:1, end:-1:1, :);
y2 = imfilter(x, h);
y2 = imfilter(y2(end:-1:1, end:-1:1, :),h);
y2 = y2(end:-1:1, end:-1:1, :);
y3 = imfilter(x, v);
y3 = imfilter(y3(end:-1:1, end:-1:1, :),v);
y3 = y3(end:-1:1, end:-1:1, :);
y = y1 + mu*(y2 + y3);
end