% 2022/10/06
% trialʱ���Ϊ�ɱ�ֵ

function disp_info=one_step(fb_offon,list,stilist_,para)

global  allimg img_rect w time0   xc yc  mountainindex cityindex endmask_texid
global  imgR sti_rect sti_R response_type time3  stitime hz

reswindow=para.reswindow;
fb_flip_time=round(time3*hz); %�׶�2��ķ���ʱ��
change_time=para.change_time;
frames_pertrial=change_time*60;

thislist=[list.imgid,0];
frame_num=[list.frame_num,frames_pertrial];

if stilist_==0
else
    stilist=stilist_;
    stilist=[stilist;0];%�������һ��endmask
end

last_texid=allimg(strcmp({allimg.type},'init')).texid;
%%
%����ÿһ֡��Ҫ���ֵ���Ϣ
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
    
    fliptime=frame_num(trial);         %ÿ��trial��֡��
    fadein=linspace(0,1,fliptime);
    fadeout=linspace(1,0,fliptime);
    sti_frame=zeros(fliptime,1);
    sti_frame(1:stitime*hz)=1;
    
    temp1=find(fadein>reswindow(2));%����ÿ��trial��Ӧ������
    temp2=find(fadein>reswindow(1));
    resframe_pertrial=zeros(1,fliptime);
    resframe_pertrial([1:temp1(1),temp2(1):fliptime])=1;
    
    if trial==1%��һ�������һ���Դη�Ӧ�����е���
        resframe_pertrial(1:temp1(1))=0;
    end
    if trial==length(thislist)
        resframe(end-length(temp2)+1:end)=0;
    end
    
    
    for flip=1:fliptime
        frame=frame+1;
        %֡ID����һ֡ID������ٷֱȣ������ٷֱ�,��Ӧ������
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
temp=find(resframe);%��Ӧ����Ӧ֡������
temp2_=[1,find(diff(temp)~=1)+1];
temp3_=[temp2_(2:end)-1,length(temp)];
for i=1:length(thislist)-1
    restrial(temp(temp2_(i):temp3_(i)))=i;
end
disp_info(:,6)=restrial';
disp_info(find(disp_info(:,6)),7)=thislist(disp_info(find(disp_info(:,6)),6));
disp_info(find([disp_info(:,3)==0]),11)=thislist;

%�þ����ڿ�ͷ��ĩβ����������initmask��endmask
%������ͼƬ�Դ���ɾ����������ͼƬ�ķ�Ӧ��

%��һ�У��Ȼ��Ƶ�ͼƬ���
%�ڶ��У�����Ƶ�ͼƬ���
%�����У��Ȼ��Ƶ�ͼƬ��͸����
%�����У�����Ƶ�ͼƬ��͸����
%�����У�0��ʾδ����Ӧ���ڣ�1��ʾ���ﷴӦ����
%�����У����ֱ�ʾ��ǰ���ڵڼ���ͼƬ�Դ�
%��7�У� ��ʾ��ǰͼƬ�Դε�ͼƬ����ֵ  Ҳ��alllist����Ϣ�����ִ�С��ʾ��ͼƬ�������ɽ���߳���

%��8�У��ڳ���ʱ���¼��ÿ�γ��ֵ�ʱ��
%��9�У� ��Ӧ���ڵĳ��ΰ�����1��ʾ����
%��10�У� ԭʼ������Ϣ

%��11��  ����ͼƬ���γ��ֵĵ�һ֡��͸����Ϊ0������Ǹ�ͼƬ������ֵ
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
        if disp_info(flip,5)-disp_info(flip-1,5)==1 %�����������ӣ��Դ�ͼƬ���� ���뷴Ӧ������
            nowtrial=thisdisp(6);   %��ǰ�ǶԵڼ����Դ�ͼƬ�ķ�Ӧ
            subresed=0; %��Ϊδ��Ӧ
            nowpic=thisdisp(7); %��ǰͼƬ������ֵ
            if nowpic>10 %��ʾ��ǰ���ִ�ɽ
                rightres=~response_type;
            elseif nowpic<=10 %��ʾ��ǰ���ֳ���
                rightres=response_type;
            end
        end
        
        if disp_info(flip,5)-disp_info(flip-1,5)==-1 %�������ּ��٣���Ӧ���ڽ��� ׼������
            if rightres==1 &&  subresed==1 %����ǰӦ�÷�Ӧ�ҷ�Ӧ�ˣ� ���¼201����ʾ����
                sendtrigger(201);
            elseif rightres==1 &&  subresed==0 %����ǰ�÷�Ӧ��û��Ӧ�� ���¼204����ʾ©��
                sendtrigger(204);
            elseif rightres==0 &&  subresed==0  %����ǰ���÷�Ӧ��û��Ӧ�����¼202����ȷ�ܳ�
                sendtrigger(202);
            elseif rightres==0 &&  subresed==1 %����ǰ���÷�Ӧ�ҷ�Ӧ�ˣ� ���¼203����ʾ�鱨
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
    
    if thisdisp(3)==0 %��ʾ��ͼƬ����  ׼�����
        if thisdisp(11)>10 %��ʾ��ǰ���ִ�ɽ�����20
            sendtrigger(20);
        elseif thisdisp(11)<=10 %��ʾ��ǰ���ֳ��У� ���10
            sendtrigger(10);
        end
    end
    
    disp_info(flip,8)=GetSecs-dispstart;
    
    [~,~,kc]=KbCheck;
    if kc(para.response_key)
        disp_info(flip,10)=1;
    end
    if thisdisp(5)==1 && subresed==0 %��ⰴ��
        if kc(para.response_key)
            disp_info(flip,9)=1;
            if rightres==1
                judge=1;
                sendtrigger(101);
            elseif rightres==0
                judge=0;
                sendtrigger(100);
            end
            fb_flip_count=fb_offon; %��ֵ��ԶΪ0 ��������ַ���
            subresed=1;
        end
    end
    
    if fb_flip_count>0
        if judge==1
            drawcenteredtext_dot(w,'��ȷ',xc,yc+imgR+100,[50 255 50],25);
        elseif judge==0
            drawcenteredtext_dot(w,'����',xc,yc+imgR+100,[255 50 50],25);
        end
        fb_flip_count=fb_flip_count+1;
    end
    if fb_flip_count>fb_flip_time
        fb_flip_count=0;
    end
    checkend;
end

%��Ӧ˲��Ľ���ж�mark��101��ȷ  100����
%��Ӧ������˲��Ľ������mark��
%����ǰӦ�÷�Ӧ�ҷ�Ӧ�ˣ� ���¼201����ʾ����
%����ǰ���÷�Ӧ��û��Ӧ�����¼202����ȷ�ܳ�
%����ǰ���÷�Ӧ�ҷ�Ӧ�ˣ� ���¼203����ʾ�鱨
%����ǰ�÷�Ӧ��û��Ӧ�� ���¼204����ʾ©��