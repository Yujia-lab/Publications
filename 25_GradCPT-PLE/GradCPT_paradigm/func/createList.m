%创建试次序列
function list = createList(trialpara)
global hz

num=trialpara.num;%试次数
targetpe_percentage=trialpara.targetpe_percentage;%目标出现概率
rand_duration=trialpara.rand_duration;%trial持续时间是否随机 1是0否
duration=trialpara.duration;%trial时长
rand_range=trialpara.rand_range;%随机范围
targetpe_rand=trialpara.targetpe_rand;%目标是否随机出现 1是0否

frame_num=duration*hz;
if frame_num~=round(frame_num)
    error('trial时长乘以刷新率必须是整数')
end

list.imgid=ceil(rand(1,num)*10);
list.imgid=getalllist(list.imgid);%相同的图片不能连续出现

%目标是否随机出现 1是0否
if targetpe_rand==1
    targetid=randperm(num,round(num*targetpe_percentage));
elseif targetpe_rand==0
    targetid=[round(1/targetpe_percentage):round(1/targetpe_percentage):num];
else
    error('目标是否随机出现 1是0否,请输入正确参数')
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