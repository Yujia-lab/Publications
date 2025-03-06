% 训练阶段
function train(para,structure_time)

global    w hz
disp_intro_mri(w,'intro1.png');
jg_start = GetSecs();

thisstep=1;
while 1
    if thisstep==1
        thislist=getalllist([randperm(10,10),randperm(10,10)+10]);
        step1(thislist,para);
        thisstep=2;
    end
    
    if thisstep==2
        rep=disp_intro(w,'intro2.png',1);
        parae2 = para;
        parae2.change_time=1.5;
        parae2.reswindow=[0.5,0.2];
        
        trialpara.num=30;%试次数
        trialpara.targetpe_percentage=0.33;%目标出现概率
        trialpara.targetpe_rand=1;%目标是否随机出现 1是0否
        trialpara.rand_duration=0;%trial持续时间是否随机 1是0否
        trialpara.duration=parae2.change_time;
        trialpara.rand_range=0;%随机范围
        list2=createList(trialpara);
        if rep==1
            thisstep=1;
            continue
        end
        
        one_step(1,list2,0,parae2);
        thisstep=3;
        disp_intro(w,'intro3.png',1);
    end
    
    if thisstep==3
        
        trialpara.num=40;%试次数
        trialpara.targetpe_percentage=0.2;%目标出现概率
        trialpara.targetpe_rand=1;%目标是否随机出现 1是0否
        trialpara.rand_duration=0;%trial持续时间是否随机 1是0否
        trialpara.duration=para.change_time;
        trialpara.rand_range=0;%随机范围
        
        list3=createList(trialpara);
        
        prac4.disp_info=one_step(0,list3,0,para);
        fliptime=round(para.change_time*hz);
        prac4.result=resultana(prac4,fliptime,0);
        prac4.cr=sum([prac4.result.judge]==1)/length([prac4.result.judge]);
        if prac4.cr>=0.9
            break
        else
            thisstep=3;
            disp_intro(w,'intro3back.png',1);
        end
    end
end

array=imread('jiegou.png');
texid=Screen('MakeTexture',w,array);
Screen('DrawTexture',w,texid);
Screen('Flip',w);
while (GetSecs-jg_start)<=structure_time
    [~,~,kc]=KbCheck;
    if kc(27)
        break;
    end
end