%�����Դ�����
function list = createList(trialpara)
global hz

num=trialpara.num;%�Դ���
targetpe_percentage=trialpara.targetpe_percentage;%Ŀ����ָ���
rand_duration=trialpara.rand_duration;%trial����ʱ���Ƿ���� 1��0��
duration=trialpara.duration;%trialʱ��
rand_range=trialpara.rand_range;%�����Χ
targetpe_rand=trialpara.targetpe_rand;%Ŀ���Ƿ�������� 1��0��

frame_num=duration*hz;
if frame_num~=round(frame_num)
    error('trialʱ������ˢ���ʱ���������')
end

list.imgid=ceil(rand(1,num)*10);
list.imgid=getalllist(list.imgid);%��ͬ��ͼƬ������������

%Ŀ���Ƿ�������� 1��0��
if targetpe_rand==1
    targetid=randperm(num,round(num*targetpe_percentage));
elseif targetpe_rand==0
    targetid=[round(1/targetpe_percentage):round(1/targetpe_percentage):num];
else
    error('Ŀ���Ƿ�������� 1��0��,��������ȷ����')
end

list.imgid(targetid)=list.imgid(targetid)+10;
list.frame_num=repmat(frame_num,1,num);

if rand_duration
    
    while 1
        rand_frame=round(rand(1,num)*rand_range*hz*2)-rand_range*hz;
        rand_frame_fix=sum(rand_frame);
        rand_frame_fix_num=abs(rand_frame_fix);
        if rand_frame_fix_num<(num/2)
            break
        end
    end
    
    rand_frame_extr_id=find((rand_frame~=max(rand_frame))&(rand_frame~=min(rand_frame)));
    rand_frame_fix_id=rand_frame_extr_id(randperm(numel(rand_frame_extr_id),rand_frame_fix_num));
    
    if rand_frame_fix > 0
        rand_frame(rand_frame_fix_id)= rand_frame(rand_frame_fix_id)-1;
    else
        rand_frame(rand_frame_fix_id)= rand_frame(rand_frame_fix_id)+1;
    end
    
    list.frame_num=list.frame_num+rand_frame;
end
end