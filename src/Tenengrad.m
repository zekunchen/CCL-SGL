function [SS] = Tenengrad(Image)
    Sx = fspecial('sobel');
    Gx = imfilter(double(Image), Sx, 'replicate', 'conv');
    Gy = imfilter(double(Image), Sx', 'replicate', 'conv');
%     [Gx, Gy] = imgradientxy(Image);
    S = sqrt(Gx.^2 + Gy.^2);
%     SS = mean2(S);
    mask=Image;

    BTG=S.*mask;
    RTG=BTG./(Image+0.00001);
    SS = mean2(RTG);


end