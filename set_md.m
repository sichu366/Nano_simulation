function set_md(T,a,time2test)
tic;[E,r,NN,NL]=md(3,4,[4,4,4],a,[1,1,1],1000,1000,20,10,5,T,time2test);toc;

h = figure;
plot3(r(:,1),r(:,2),r(:,3),'o')
grid on;

str = 'Final position,T=';
Tstr = num2str(T);
str = strcat(str,Tstr);
str = strcat(str,',Box side len=');
boxStr = num2str(a);
str = strcat(str,boxStr);
title(str,'fontsize',18);
basepath = pwd;
path_save = strcat(basepath,num2str(time2test));
path_save = strcat(path_save,str);
path_save = strcat(path_save,'.jpg');
saveas(h,path_save);

h = figure;
plot(NN,'k-o')
xlim([0 260])
grid on;
str = 'Num of Neighbors,T=';
Tstr = num2str(T);
str = strcat(str,Tstr);
str = strcat(str,',Box side len=');
boxStr = num2str(a);
str = strcat(str,boxStr);

title(str,'fontsize',18);
path_save = strcat(basepath,num2str(time2test));
path_save = strcat(path_save,str);
path_save = strcat(path_save,'.jpg');
saveas(h,path_save);


E_mean=mean(E);
E_std=std(E);
t=(1:50)*0.1; % ps
h = figure;
plot(t,E(:,2),'^-',t,E(:,1),'v-',t,E(:,3),'o-')
xlabel('Time (ps)','fontsize',18);
ylabel('Energy (eV)','fontsize',18);
legend({'Kinetic','Potential','Total'},'position',[0.4,0.6,0.25,0.1]); 
set(gca,'fontsize',18);
legend('boxoff') 

str = 'Energy (eV),T=';
Tstr = num2str(T);
str = strcat(str,Tstr);
str = strcat(str,',Box side len=');
boxStr = num2str(a);
str = strcat(str,boxStr);

path_save = strcat(basepath,num2str(time2test));
path_save = strcat(path_save,str);
path_save = strcat(path_save,'.jpg');
saveas(h,path_save);