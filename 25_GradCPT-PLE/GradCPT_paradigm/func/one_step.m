% 2022/10/06
% trial时间改为可变值

function disp_info=one_step(fb_offon,list,stilist_,para)

global  allimg img_rect w time0   xc yc  mountainindex cityindex endmask_texid
global  imgR sti_rect sti_R response_type time3  stitime hz

reswindow=para.reswindow;
fb_flip_time=round(time3*hz); %阶段2里的反馈时长
change_time=para.change_time;
frames_pertrial=change_time*60;

thislist=[list.imgid,0];
frame_num=[list.frame_num,frames_pertrial];

if stilist_==0
else
    stilist=stilist_;
    stilist=[stilist;0];%最后增加一个endmask
end

last_texid=allimg(strcmp({allimg.type},'init')).texid;
%%
%计算每一帧需要呈现的信息
disp_info=[];
frame=0;
for trial=1:length(thislist)
    thistrial=thislist(trial);
    if thistrial>10
        this_texid=allimg(mountainindex(thistrial-10)).texid;
    elseif thistrial==0
        this_texid=endmask_texid;
    elseif thistrial<=10
        this_texid=allimg(cityindex(thistrial)).texid;
    end
    
    fliptime=frame_num(trial);         %每个trial的帧数
    fadein=linspace(0,1,fliptime);
    fadeout=linspace(1,0,fliptime);
    sti_frame=zeros(fliptime,1);
    sti_frame(1:stitime*hz)=1;
    
    temp1=find(fadein>reswindow(2));%生成每个trial反应窗序列
    temp2=find(fadein>reswindow(1));
    resframe_pertrial=zeros(1,fliptime);
    resframe_pertrial([1:temp1(1),temp2(1):fliptime])=1;
    
    if trial==1%第一个和最后一个试次反应窗序列调整
        resframe_pertrial(1:temp1(1))=0;
    end
    if trial==length(thislist)
        resframe(end-length(temp2)+1:end)=0;
    end
    
    
    for flip=1:fliptime
        frame=frame+1;
        %帧ID，下一帧ID，淡入百分比，淡出百分比,反应窗序列
        disp_info(frame,:)=[this_texid,last_texid,fadein(flip),fadeout(flip),resframe_pertrial(flip)];
        resframe(frame)=resframe_pertrial(flip);
        
        if stilist_==0
        else
            sti_frame_list(frame)=sti_frame(flip)*stilist(trial);
        end
    end
    
    last_texid=this_texid;
end
%%

restrial=resframe;
temp=find(resframe);%反应窗对应帧的索引
temp2_=[1,find(diff(temp)~=1)+1];
temp3_=[temp2_(2:end)-1,length(temp)];
for i=1:length(thislist)-1
    restrial(temp(temp2_(i):temp3_(i)))=i;
end
disp_info(:,6)=restrial';
disp_info(find(disp_info(:,6)),7)=thislist(disp_info(find(disp_info(:,6)),6));
disp_info(find([disp_info(:,3)==0]),11)=thislist;

%该矩阵在开头和末尾处，加入了initmask和endmask
%并且在图片试次里删掉了这两个图片的反应窗

%第一列，先绘制的图片句柄
%第二列，后绘制的图片句柄
%第三列，先绘制的图片的透明度
%第四列，后绘制的图片的透明度
%第五列，0表示未到反应窗口，1表示到达反应窗口
%第六列，数字表示当前处于第几个图片试次
%第7列， 表示当前图片试次的图片索引值  也即alllist的信息，数字大小表示了图片的类别是山或者城市

%第8列，在呈现时候记录了每次呈现的时刻
%第9列， 反应窗内的初次按键，1表示按键
%第10列， 原始按键信息

%第11列  当新图片初次出现的第一帧（透明度为0），标记该图片的索引值
%%
last_texid=allimg(strcmp({allimg.type},'init')).texid;
Screen('DrawTexture',w,last_texid,[],img_rect);
if stilist_==0
else
    Screen('Frameoval',w,0,sti_rect,sti_R);
end
Screen('Flip',w);
WaitSecs(time0);

dispstart=GetSecs;
subresed=0;
fb_flip_count=0;
for flip=1:size(disp_info,1)
    thisdisp=disp_info(flip,:);
    
    if flip>1
        if disp_info(flip,5)-disp_info(flip-1,5)==1 %发生数字增加，试次图片更新 进入反应窗口内
            nowtrial=thisdisp(6);   %当前是对第几个试次图片的反应
            subresed=0; %记为未反应
            nowpic=thisdisp(7); %当前图片的索引值
            if nowpic>10 %表示当前呈现大山
                rightres=~response_type;
            elseif nowpic<=10 %表示当前呈现城市
                rightres=response_type;
            end
        end
        
        if disp_info(flip,5)-disp_info(flip-1,5)==-1 %发生数字减少，反应窗口结束 准备结算
            if rightres==1 &&  subresed==1 %若当前应该反应且反应了， 则记录201，表示击中
                sendtrigger(201);
            elseif rightres==1 &&  subresed==0 %若当前该反应且没反应， 则记录204，表示漏报
                sendtrigger(204);
            elseif rightres==0 &&  subresed==0  %若当前不该反应且没反应，则记录202，正确拒斥
                sendtrigger(202);
            elseif rightres==0 &&  subresed==1 %若当前不该反应且反应了， 则记录203，表示虚报
                sendtrigger(203);
            end
        end
    end
    
    Screen('DrawTexture',w,thisdisp(1),[],img_rect,[],[],thisdisp(3));
    Screen('DrawTexture',w,thisdisp(2),[],img_rect,[],[],thisdisp(4));
    if stilist_==0
    else
        Screen('Frameoval',w,sti_frame_list(flip)*255,sti_rect,sti_R);
    end
    Screen('Flip',w);
    
    if thisdisp(3)==0 %表示新图片出现  准备打标
        if thisdisp(11)>10 %表示当前呈现大山，打标20
            sendtrigger(20);
        elseif thisdisp(11)<=10 %表示当前呈现城市， 打标10
            sendtrigger(10);
        end
    end
    
    disp_info(flip,8)=GetSecs-dispstart;
    
    [~,~,kc]=KbCheck;
    if kc(para.response_key)
        disp_info(flip,10)=1;
    end
    if thisdisp(5)==1 && subresed==0 %检测按键
        if kc(para.response_key)
            disp_info(flip,9)=1;
            if rightres==1
                judge=1;
                sendtrigger(101);
            elseif rightres==0
                judge=0;
                sendtrigger(100);
            end
            fb_flip_count=fb_offon; %该值永远为0 ，不会呈现反馈
            subresed=1;
        end
    end
    
    if fb_flip_count>0
        if judge==1
            drawcenteredtext_dot(w,'正确',xc,yc+imgR+100,[50 255 50],25);
        elseif judge==0
            drawcenteredtext_dot(w,'错误',xc,yc+imgR+100,[255 50 50],25);
        end
        fb_flip_count=fb_flip_count+1;
    end
    if fb_flip_count>fb_flip_time
        fb_flip_count=0;
    end
    checkend;
end

%反应瞬间的结果判断mark，101正确  100错误
%反应窗结束瞬间的结果结算mark：
%若当前应该反应且反应了， 则记录201，表示击中
%若当前不该反应且没反应，则记录202，正确拒斥
%若当前不该反应且反应了， 则记录203，表示虚报
%若当前该反应且没反应， 则记录204，表示漏报