function  step1(thislist,para)

global  allimg img_rect w time0 time1 time2 xc yc  mountainindex cityindex

last_texid=allimg(strcmp({allimg.type},'init')).texid;
Screen('DrawTexture',w,last_texid,[],img_rect);
Screen('Flip',w);
WaitSecs(time0);

trial=1;
while 1
    thistrial=thislist(trial);
    if thistrial>10
        this_texid=allimg(mountainindex(thistrial-10)).texid;
        this_type=allimg(mountainindex(thistrial-10)).type;
        this_name=allimg(mountainindex(thistrial-10)).name;
    else
        this_texid=allimg(cityindex(thistrial)).texid;
        this_type=allimg(cityindex(thistrial)).type;
        this_name=allimg(cityindex(thistrial)).name;
    end
    Screen('DrawTexture',w,this_texid,[],img_rect);
    Screen('Flip',w);
    subresed=0;
    resstart=GetSecs;
    while 1
        if GetSecs-resstart>time1
            break;
        end
        [~,~,kc]=KbCheck;
        if kc(para.response_key)
            subresed=1;
            break;
        end
        checkend;
    end
    if subresed==1 && thistrial<=10  %����ͼƬ  ����
        judge=1;
        drawcenteredtext_dot(w,'������',xc,yc,[50,255,50],30);
    elseif subresed==0 && thistrial>10   %��ɽͼƬ û����
        drawcenteredtext_dot(w,'������',xc,yc,[50,255,50],30);
        judge=1;
    elseif subresed==0 && thistrial<10   %����ͼƬ û����
        drawcenteredtext_dot(w,'�����ˣ���������Ҫ����',xc,yc,[250,55,50],30);
        judge=0;
    elseif subresed==1 && thistrial>=10   %����ͼƬ û����
        drawcenteredtext_dot(w,'�����ˣ�������ɽ���ܰ���',xc,yc,[250,50,50],30);
        judge=0;
    end
    Screen('Flip',w);
    WaitSecs(time2);    
    if judge==1
        trial=trial+1;
    end
    if  trial>length(thislist)
        break;
    end
end