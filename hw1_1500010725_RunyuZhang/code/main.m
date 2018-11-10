function main (image_path , lambda_weight , ...
kernel_size , gaussian_sigma , result_image_path, opts )
% ---- Inputs:
% image_path
% lambda_weight
% kernel_size
% gaussian_sigma
% result_image_path
% opts: A struct for optional inputs, including:
%   mu: The coefficent for the second order penalizing term, default 1
%   tmax: The maximum allowable iteration time, default 30
%   type: The boundary condition for the blur kernel, has choice 'cyclic' &
%   'noncyclic', default 'cyclic'
if nargin <1; image_path='image/kamiya.jpg' ; end % ����ͼƬ��
if nargin <2; lambda_weight =1; end % lambda
if nargin <3; kernel_size =15; end % ����˴�С
if nargin <4; gaussian_sigma = 1.5 ; end % ��˹�˷���
if nargin <5; result_image_path = ...
'image/kamiya_result.jpg' ; end % ���ͼƬ��
if nargin <6; mu = 1; tmax = 30; type = 'cyclic';
else
    if isfield(opts,'mu'); mu = opts.mu; else mu = 1;end
    if isfield(opts,'tmax'); tmax = opts.tmax; else tmax = 30;end
    if isfield(opts, 'type'); type = opts.type; else type = 'cyclic';end
end
% your code , e . g .
kernel = fspecial('gaussian',[kernel_size, kernel_size] , gaussian_sigma);
u = imread (image_path); % ��ȡͼƬ
u = im2double(u);
switch type
    case 'cyclic'
        mysolver = @tv_deblur_cyclic;
    case 'noncyclic'
        mysolver = @tv_deblur_noncyclic;
end
u = mysolver ( u , kernel , lambda_weight, mu, tmax); % �Լ�ʵ���������
h = figure ;
imshow(u)
print( h , result_image_path ,'-dpng') % ����洢������
end