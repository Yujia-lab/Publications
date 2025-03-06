function allimg=readallimg

global w

allimg=dir('pic/*.png');
for i=1:length(allimg)
    thisname=allimg(i).name;
    if strfind(thisname,'city')
        allimg(i).type='city';
    elseif strfind(thisname,'mountain')
        allimg(i).type='mountain';
    else
        allimg(i).type='init';
    end
    array=imread(['pic/',thisname]);
    [img_h,img_w,channel]=size(array);
    if img_h>=img_w
        array=array(round((img_h-img_w)/2)+1:round((img_h-img_w)/2)+img_w,:,:);
    else
        array=array(:,round((img_w-img_h)/2) +1: round((img_w-img_h)/2)+img_h,:);
    end
    [img_s,~,~]=size(array);
    lintemp=linspace(-img_s/2,img_s/2,img_s);
    [X,Y]=meshgrid(lintemp,lintemp);
    mask=sqrt(X.^2+Y.^2);
    
    R=mask(1,round(size(mask,1)/2));
    index1=mask<=R;
    index2=mask>R;
    mask(index1)=255;
    mask(index2)=0;
    array(:,:,size(array,3)+1)=mask;
    allimg(i).array=array;
    allimg(i).texid=Screen('MakeTexture',w,array);
end
