function disp_intro_mri(w,array_name)
%接收S信号开始实验
array=imread(array_name);
texid=Screen('MakeTexture',w,array);
Screen('DrawTexture',w,texid);
Screen('Flip',w);
while 1    
    [~,kc]=KbStrokeWait();
    if kc(83)
        break;
    end    
end
Screen('Close',texid);
