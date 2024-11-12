function [degradedIm, H_full] = addblur(h, img)
    [Ny,Nx,Nz]=size(img);%获取输入图像 img 的大小，Ny 是图像的高度，Nx 是图像的宽度，Nz 是图像的通道数（如果是彩色图像，则为3
    H_full = zeros(Ny,Nx);%创建一个大小为图像大小的全零矩阵 H_full，用于存储卷积核
    H_full(1:size(h,1), 1:size(h,2)) = h;%将输入的卷积核 h 复制到 H_full 中。如果 h 的大小小于 H_full，则 h 将被放置在左上角
    center=[(size(h,1)+1)/2,(size(h,2)+1)/2];%计算卷积核 h 的中心坐标，用于后续的频域移位
    blur_A=padPSF(H_full,[Ny,Nx]);
    blur_matrix_trans=fft2(circshift(blur_A,1-center));
    degradedIm = real(ifft2(fft2(img).*blur_matrix_trans));
    H_full = circshift(blur_A,1-center);
end