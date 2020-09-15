%% Porkchop Plot

clear C3

tof_min=365; %days
ToF=86400*tof_min; %seconds
d_dep=[1:15:601]; %days
d_arr=linspace(d_dep(1)+4*tof_min,d_dep(end)+8*tof_min,length(d_dep)); %days
m=1;
for i=d_dep
    n=1;
    d_arr=round(d_arr);
    for j=d_arr
        tof=j-i;
        ToF=86400*tof;
        r1=system_N(i,19:21);
        r2=system_N(tof,37:39);
        [v1,v2] = Lambert_solver(r1,r2,ToF,'long',mu(1));
        dv1=norm(system_N(i,22:24)-v1);
        dv2=norm(system_N(tof,40:42)-v2);
        C3(m,n)=dv1^2;
        n=n+1;
    end
    m=m+1;
end

[Y,I] = min(C3)
figure('name','Porkchop Plot')
contour(d_dep,d_arr,C3,'ShowText','on')
hold on
title('Earth to Saturn Porkchop Plot')
% for k=1:length(I)
% plot(d_dep(k),d_arr(I(k)),'o')
% end
plot(d_dep(end),d_arr(I(end)),'o')
xlabel('Depature Date - days after 01/01/2020')
ylabel('Arrival Date - days after 01/01/2020')
text(d_dep(end),d_arr(I(end)),strcat('\leftarrow','minimum C3 =',num2str(Y(end))),'HorizontalAlignment','left')
legend('C3')