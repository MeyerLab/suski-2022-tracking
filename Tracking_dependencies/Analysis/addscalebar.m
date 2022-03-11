function RGB=addscalebar(RGB,startpos,pixdims)
    h=size(RGB,1);
    w=size(RGB,2);
    barimg=zeros(h,w);
    barimg(startpos(1):startpos(1)+pixdims(1),startpos(2):startpos(2)+pixdims(2))=1;
    RGB(:,:,1)=RGB(:,:,1)+barimg;
    RGB(:,:,2)=RGB(:,:,2)+barimg;
    RGB(:,:,3)=RGB(:,:,3)+barimg;
    RGB(RGB>1)=1;
end