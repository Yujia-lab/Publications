function rep=disp_intro(w,array_name,type)
 %typeΪ1����ʾֻ�ܰ��س���
 %typeΪ2����ʾ���԰��س���R��

 
array=imread(array_name);
texid=Screen('MakeTexture',w,array);
while 1
    Screen('DrawTexture',w,texid);
    Screen('Flip',w);
    [~,~,kc]=KbCheck;
    
    if kc(82) && type==2
        rep=1;
        break;
    end
    
    if kc(50)% 2@
        rep=0;
        break;
    end
    
end
Screen('Close',texid);
KbReleaseWait();