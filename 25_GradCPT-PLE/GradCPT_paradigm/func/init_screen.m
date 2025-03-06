function [w,wrect,hz,xc,yc]=init_screen(bgc,skip)

Screen('CloseAll');
Screen('Preference','SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', skip); 
screens=Screen('Screens');  screen_num=max(screens);  
[w,wrect]=Screen('OpenWindow',screen_num,bgc,[]);
ListenChar(2) ;  %�رձ༭���Լ��̼���
HideCursor;
hz=FrameRate(w);
hz=round(hz);
[xc,yc]=WindowCenter(w);
AssertOpenGL;
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Priority(MaxPriority(w));
KbName('UnifyKeyNames');
Screen('TextSize',w,25);
Screen('TextFont',w,'-:lang=zh-cn');
rng('Shuffle');