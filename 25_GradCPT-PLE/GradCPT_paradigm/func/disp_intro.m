function rep=disp_intro(w,array_name,type)
 %type为1，表示只能按回车键
 %type为2，表示可以按回车和R键

 
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