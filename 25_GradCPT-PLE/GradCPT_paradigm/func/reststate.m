function reststate(rest_time)

global    w    xc yc   imgR bgc

    disp_intro_mri(w,'jingxi.png');
    Screen('FillRect',w,bgc);
    Screen('DrawLines',w,[xc,xc,xc+0.7*imgR,xc-0.7*imgR;yc+0.7*imgR,yc-0.7*imgR,yc,yc],7,0);
    Screen('Flip',w);
    jx_start = GetSecs;
    while (GetSecs-jx_start)<=rest_time
        [~,~,kc]=KbCheck;
        if kc(27)
            break;
        end
    end
end