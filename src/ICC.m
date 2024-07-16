% Image correlation Coeffcient ICC
function out=ICC(sig,sgt)
    up=sum((sig-mean(sig(:))).*(sgt-mean(sgt(:))));
    do=norm(sig-mean(sig(:)))*norm(sgt-mean(sgt(:)));
    out=up/do;
end