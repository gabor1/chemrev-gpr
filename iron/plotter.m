%%%%

dr_model

Stress_model = K*tau_P*eVtoGPa;
DH_model = Un;

tauqp = 0:0.01:1;
%DHqp = (1-tauqp.^(3/4)).^(5/4);
%DHqp = (1-tauqp.^(0.75)).^(1.25);
DHqp = (1-tauqp.^(0.77)).^(1.30);

rect=[300 300 1200 1000];

p=[];
hN=[];

figure(1000); clf

h1 = subplot(1,2,1)
set(gcf,'OuterPosition',rect);
p(1)=plot(GAP_40b_0MPa_F_32R(:,1),(GAP_40b_0MPa_F_32R(:,2)-GAP_40b_0MPa_F_32R(1,2)),'k .-','MarkerSize',30); hold on
p(2)=plot(GAP_40b_100MPa_F_32R(:,1),(GAP_40b_100MPa_F_32R(:,2)-GAP_40b_100MPa_F_32R(1,2)),'r .-','MarkerSize',30);
p(3)=plot(GAP_40b_200MPa_F_32R(:,1),(GAP_40b_200MPa_F_32R(:,2)-GAP_40b_200MPa_F_32R(1,2)),'r .--','MarkerSize',30);
p(4)=plot(GAP_40b_400MPa_F_32R(:,1),(GAP_40b_400MPa_F_32R(:,2)-GAP_40b_400MPa_F_32R(1,2)),'r .-.','MarkerSize',30);
p(5)=plot(GAP_40b_500MPa_F_32R(:,1),(GAP_40b_500MPa_F_32R(:,2)-GAP_40b_500MPa_F_32R(1,2)),'m .-','MarkerSize',30);
p(6)=plot(GAP_40b_600MPa_F_32R(:,1),(GAP_40b_600MPa_F_32R(:,2)-GAP_40b_600MPa_F_32R(1,2)),'m .--','MarkerSize',30);
p(7)=plot(GAP_40b_800MPa_F_32R(:,1),(GAP_40b_800MPa_F_32R(:,2)-GAP_40b_800MPa_F_32R(1,2)),'m .-.','MarkerSize',30);
p(8)=plot(GAP_40b_1000MPa_F_32R(:,1),(GAP_40b_1000MPa_F_32R(:,2)-GAP_40b_1000MPa_F_32R(1,2)),'c .-','MarkerSize',30);
p(9)=plot(GAP_40b_1200MPa_F_32R(:,1),(GAP_40b_1200MPa_F_32R(:,2)-GAP_40b_1200MPa_F_32R(1,2)),'c .--','MarkerSize',30);
p(10)=plot(GAP_40b_1400MPa_F_32R(:,1),(GAP_40b_1400MPa_F_32R(:,2)-GAP_40b_1400MPa_F_32R(1,2)),'c .-.','MarkerSize',30);
p(11)=plot(GAP_40b_1500MPa_F_32R(:,1),(GAP_40b_1500MPa_F_32R(:,2)-GAP_40b_1500MPa_F_32R(1,2)),'b .-','MarkerSize',30);
p(12)=plot(GAP_40b_1600MPa_F_32R(:,1),(GAP_40b_1600MPa_F_32R(:,2)-GAP_40b_1600MPa_F_32R(1,2)),'b .--','MarkerSize',30);
p(13)=plot(GAP_40b_1800MPa_F_32R(:,1),(GAP_40b_1800MPa_F_32R(:,2)-GAP_40b_1800MPa_F_32R(1,2)),'b .-.','MarkerSize',30);
p(14)=plot(GAP_40b_2000MPa_F_32R(:,1),(GAP_40b_2000MPa_F_32R(:,2)-GAP_40b_2000MPa_F_32R(1,2)),'g .-','MarkerSize',30);

text('FontSize',25,'HorizontalAlignment','center')

title('a','Fontname','Arial','fontsize',40,'fontweight','bold','Position',[0.01 1.01 0]);

%title('Enthalpy change (eV) vs reaction coordinate as fct. of applied stress','Interpreter','Latex','fontsize',22,'fontweight','bold');
y=ylabel('\Delta H(\tau) (eV)');
x=xlabel('Reaction coordinate \xi (-)');

l=legend('0 GPa','0.1 GPa','0.2 GPa','0.4 GPa','0.5 GPa','0.6 GPa','0.8 GPa','1.0 GPa','1.2 GPa','1.4 GPa','1.5 GPa','1.6 GPa','1.8 GPa','2.0 GPa','Location','NorthWest');

%set(l,'Interpreter','Latex');
set(l,'Fontname','Times','fontsize',17);
set(y,'Fontname','Times','fontsize',25);
set(x,'Fontname','Times','fontsize',25);
set(gca,'Fontname','Times','FontSize',20);

for i=1:size(p,2)
set(p(i),'LineWidth',2);
end

axis([0 0.5 0.0 1.0]);

grid on

h2 = subplot(1,2,2)
set(gcf,'OuterPosition',rect);
hN(1)=plot(Stress_model,DH_model,'b -'); hold on
hN(2)=plot(tauqp*tau_P*eVtoGPa,DHqp*Ek,'k -');
hN(3)=plot(Barriers_Long(1,1)/1000,Barriers_Long(1,2),'k .','MarkerSize',30);
hN(4)=plot(Barriers_Long(2:4,1)/1000,Barriers_Long(2:4,2),'r .','MarkerSize',30);
hN(5)=plot(Barriers_Long(5:7,1)/1000,Barriers_Long(5:7,2),'m .','MarkerSize',30);
hN(6)=plot(Barriers_Long(8:10,1)/1000,Barriers_Long(8:10,2),'c .','MarkerSize',30);
hN(7)=plot(Barriers_Long(11:13,1)/1000,Barriers_Long(11:13,2),'b .','MarkerSize',30);
hN(8)=plot(Barriers_Long(14,1)/1000,Barriers_Long(14,2),'g .','MarkerSize',30);

text('FontSize',25,'HorizontalAlignment','center')

title('b','Fontname','Arial','fontsize',40,'fontweight','bold','Position',[0.01 1.01 0]);

%title('Enthalpy change (eV) vs applied stress','Interpreter','Latex','fontsize',22,'fontweight','bold');
%y=ylabel('\Delta H(\tau) (eV)');
x=xlabel('Applied stress \tau (GPa)');

l=legend([hN(1) hN(2)],{'\DeltaH^*(\tau), LT model','\DeltaH^*(\tau), Kocks law (p=0.77, q=1.3)',},'Location','NorthEast');

%l=legend('\DeltaH^*(\tau), GAP','\DeltaH^*(\tau), LT model','\DeltaH^*(\tau), Kocks law (p=0.77, q=1.3)','Location','NorthEast');
set(l,'Fontname','Times','fontsize',17);
%set(y,'Fontname','Times','fontsize',25);
set(x,'Fontname','Times','fontsize',25);
set(gca,'Fontname','Times','FontSize',20);
set(gca,'YTickLabel',[]);

for i=1:size(hN,2)
set(hN(i),'LineWidth',2);
end

axis([0 2.1 0 1]);

grid on

p1 = get(h1, 'Position');
p2 = get(h2, 'Position');

p1(1) = 0.08;
p1(3) = 0.47;
set(h1, 'pos', p1);

p2(3) = 0.4;
set(h2, 'pos', p2);
