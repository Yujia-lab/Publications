function all_gpm = plot_age_gs(age,GS_spec,f);

%% display GS picture


frac_GS1 = [];
frac_GS2 = [];
frac_GS3 = [];
frac_GS4 = [];
frac_GS5 = [];
frac_GS6 = [];



for isub = 1:300
    if age(isub) <= 30
        frac_GS1 = [frac_GS1,GS_spec(:,isub)];
    end
    if age(isub) >= 31 & age(isub) <= 40
        frac_GS2 = [frac_GS2, GS_spec(:,isub)];
    end
    if age(isub) >= 41 & age(isub) <= 50
        frac_GS3 = [frac_GS3, GS_spec(:,isub)];
    end
    if age(isub) >= 51 & age(isub) <= 60
        frac_GS4 = [frac_GS4, GS_spec(:,isub)];
    end
    if age(isub) >= 61 & age(isub) <= 70
        frac_GS5 = [frac_GS5, GS_spec(:,isub)];
    end
    if age(isub) >= 71 
        frac_GS6 = [frac_GS6, GS_spec(:,isub)];
    end

end

gpm1 = mean(frac_GS1,2);
gpm2 = mean(frac_GS2,2);
gpm3 = mean(frac_GS3,2);
gpm4 = mean(frac_GS4,2);
gpm5 = mean(frac_GS5,2);
gpm6 = mean(frac_GS6,2);


all_gpm = [gpm1,gpm2,gpm3,gpm4,gpm5,gpm6];

plot(f,gpm1,'color',[1,1,0],'linew',2);
hold on
plot(f,gpm2,'color',[0.9,0.9,0],'linew',2);
plot(f,gpm3,'color',[0.8,0.8,0],'linew',2);
plot(f,gpm4,'color',[0.7,0.7,0],'linew',2);
plot(f,gpm5,'color',[0.65,0.65,0],'linew',2);
plot(f,gpm6,'color',[0.6,0.6,0],'linew',2);

hold off
