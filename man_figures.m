
% DO observations
if (~exist('happihavainnot', 'var'))
    if (~exist('Obs_O2', 'var'))
        Obs_O2=xlsread...
            ('I_O\KJ_Obs_O2.xls','oxygen');
        Obs_O2(1:2, :)=[ ];
        O2Obs=[datenum(Obs_O2(:,1:3)) - datenum(year,1,1),Obs_O2(:,4), Obs_O2(:,5)];
    end
    [d1,d2]=size(Obs_O2);
    maara = (d1)/23;
    happihavainnot=NaN*zeros(maara,3+23);
    paivavek = NaN*zeros(maara,6);
    paivanum = NaN*zeros(maara,1);
    
    for i=1:maara
        paivavek(i,:) = [Obs_O2(23*(i-1)+1,1),Obs_O2(23*(i-1)+1,2),Obs_O2(23*(i-1)+1,3),0,0,0];
        happihavainnot(i,1:3) = [paivavek(i,1),paivavek(i,2),paivavek(i,3)];
        for j=0:22
            happihavainnot(i,j+4)=Obs_O2(1+23*(i-1)+j,5);
        end
        
    end
end

%CO2 observations

if (~exist('CO2Obs', 'var'))
    [NumMx,StrMx] = xlsread...
        ('I_O\Kuiva_Obs_CO2.xls','CO2');
    dumda = datevec(StrMx(:,1),'dd.mm.yyyy');
    CO2depths = NumMx(1,:);
    NumMx(1,:) = [];
    Obs_CO2 = NaN*ones(size(dumda,1)*length(CO2depths),5);
    for k = 1:size(dumda,1)
        for m = 1:length(CO2depths)
            Obs_CO2((k-1)*length(CO2depths)+m,1:3) = dumda(k,1:3);
            Obs_CO2((k-1)*length(CO2depths)+m,4) = CO2depths(m);
            Obs_CO2((k-1)*length(CO2depths)+m,5) = NumMx(k,m);
        end
        
    end
    CO2Obs=[datenum(Obs_CO2(:,1:3)) - datenum(year,1,1),Obs_CO2(:,4), Obs_CO2(:,5)];
    CO2Obs_korjattu = CO2Obs;
    CO2Obs_korjattu(CO2Obs_korjattu(:,3)<0.5,3) = NaN;
end

if (~exist('CO2ObsAuto', 'var'))
    [CO2ObsAuto,~] = xlsread...
        ('I_O\Kuiva_Obs_CO2.xls','CO2auto');
end

%============================ TEMPERATURE ==================================
%%{

outpos_PPNW = [0, 100, 1300, 900];
font_big_PPNW = 16;
font_med_PPNW = 14;

vakx = 0.10;
vaky = 0.05;
lev = 0.76;
kork = 0.38;
oikealle = 0;
ylos = 0.50;

ionoff = [datenum([2013 5 2])- datenum([year 1 1]),datenum([2013 11 25])- datenum([year 1 1]),...
    datenum([2014 4 15])- datenum([year 1 1]),datenum([2014 12 23])- datenum([year 1 1])];
ionoff_obs = [datenum([2013 5 1]),datenum([2013 11 27]),datenum([2014 4 12]),datenum([2015 12 18])];

figure(220)
clf
set(gcf, 'Outerposition',[400 0 800 1000])

subplot(211)

pcolor(tt_mod,zz+dz/2,Tzt)
shading interp
hold on
hhhh = plot(ionoff,zeros(1,4),'kv','MarkerSize',8);
axis ij
datetick('x','mmm');
set(gca,'xticklabel',[]);
set(gca,'ylim',[0 13.5]);
set(gca,'xlim',[datenum([2013 1 1])-datenum([year 1 1]) datenum([2014 12 31])-datenum([year 1 1])])
ylabel('Depth (m)','Fontsize',font_med_PPNW)
caxis([0 27]);
H = colorbar('location','eastoutside');
set(H,'position',[0.91 0.295 0.04 kork])
set(gca,'fontsize',font_med_PPNW)
set(gca,'Position',[vakx vaky+ylos lev kork])
text(-60,-1.1,'Simulated','Fontsize',font_big_PPNW)
text(413,7.6,'^oC','Fontsize',font_med_PPNW)

annotation('line',[0.108,0.859],[0.53,0.53],'LineWidth',2);
annotation('line',[0.108,0.108],[0.52,0.54],'LineWidth',2);
annotation('line',[0.859,0.859],[0.52,0.54],'LineWidth',2);
annotation('line',[0.480,0.480],[0.52,0.54],'LineWidth',2);
text(datenum([2013 7 1])-datenum([year 1 1]),15.0,'Calibration period','HorizontalAlignment','center','FontSize',font_med_PPNW)
text(datenum([2014 7 1])-datenum([year 1 1]),15.0,'Validation period','HorizontalAlignment','center','FontSize',font_med_PPNW)

subplot(212)

daygrid = datenum([2013 1 1]):datenum([2014 12 31]);
pcolor(daygrid,depths',numMx')
shading interp
hold on
plot(ionoff_obs,zeros(1,4),'kv','MarkerSize',8);
axis ij
caxis([0 27]);
datetick('x','mmm')
set(gca,'ylim',[0 13.5]);
set(gca,'xlim',[datenum([2013 1 1]) datenum([2014 12 31])])
ylabel('Depth (m)','Fontsize',font_med_PPNW)
set(gca,'Position',[vakx vaky lev kork])
text(datenum([2013 11 1]),-1.0,'Observed','Fontsize',font_big_PPNW)
set(gca,'fontsize',font_med_PPNW)

set(gcf,'Color','w')
%}

%=========================== DO ===================================

outpos = [100, 100, 1100, 900]; %1283, 978
tlims_sim = [datenum([2013 1 1]) datenum([2015 1 1])];

happitasot = [1 5 9];
ht2 = happitasot;

O2zt_mod=interp1(zz+dz/2,O2zt,ht2);

figure(420) %Mallinnetut hapet
clf
kertoja = 1000/32;
set(gcf, 'Outerposition',outpos)
subplot(311)

inx=find(O2Obs(:,2)==happitasot(1));
plot(O2Obs(inx,1)+datenum(year,1,1),kertoja*O2Obs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),kertoja*O2zt_mod(1,:)./1000,'-b','LineWidth',1.5);
set(gca,'ylim',kertoja*[0 15]);
%ylabel('mg l^{-1}','fontsize',20)
datetick('x','mmm')
%grid on
set(gca,'fontsize',font_big_PPNW)
set(gca,'xticklabel',[]);
set(gca,'xlim',tlims_sim)
klpois = plot([datenum([2013 1 8]),datenum([2014 12 31])],kertoja*[15.8,15.8],'k','LineWidth',1.5);
klpois2 = plot([datenum([2013 1 8]),datenum([2013 1 8])],kertoja*[15.8,17],'k','LineWidth',1.5);
klpois3 = plot([datenum([2014 1 1]),datenum([2014 1 1])],kertoja*[15.8,17],'k','LineWidth',1.5);
klpois4 = plot([datenum([2014 12 31]),datenum([2014 12 31])],kertoja*[15.8,17],'k','LineWidth',1.5);
set(klpois,'Clipping','off')
set(klpois2,'Clipping','off')
set(klpois3,'Clipping','off')
set(klpois4,'Clipping','off')
hold off
text(datenum([2013 7 1]),kertoja*17.5,'Calibration period','HorizontalAlignment','center','FontSize',font_med_PPNW)
text(datenum([2014 7 1]),kertoja*17.5,'Validation period','HorizontalAlignment','center','FontSize',font_med_PPNW)
text(735250,kertoja*13.5,'(a) 1 m','Fontsize',font_med_PPNW)
set(gca,'Position',[0.09 0.67 0.880 0.2600])

subplot(312)

inx=find(O2Obs(:,2)==happitasot(2));
plot(O2Obs(inx,1)+datenum(year,1,1),kertoja*O2Obs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),kertoja*O2zt_mod(2,:)./1000,'-b','LineWidth',1.5);
set(gca,'ylim',kertoja*[0 15]);
ylabel('DO (mmol m^{-3})','fontsize',font_big_PPNW)
datetick('x','mmm')
%grid on
set(gca,'fontsize',font_big_PPNW)
set(gca,'xticklabel',[]);
set(gca,'xlim',tlims_sim)
text(735250,kertoja*13.5,'(b) 5 m','Fontsize',font_med_PPNW)
set(gca,'Position',[0.09 0.36 0.880 0.2600])

subplot(313)

inx=find(O2Obs(:,2)==happitasot(3));
plot(O2Obs(inx,1)+datenum(year,1,1),kertoja*O2Obs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),kertoja*O2zt_mod(3,:)./1000,'-b','LineWidth',1.5);
set(gca,'ylim',kertoja*[0 15]);

datetick('x','mmm')
set(gca,'fontsize',font_big_PPNW)
set(gca,'xlim',tlims_sim)
set(gca,'Position',[0.09 0.05 0.880 0.2600])
HO2 = legend('Measurement','Simulation');
set(HO2,'Location','SouthOutside','Orientation','horizontal','Box','on');
set(HO2,'fontsize',font_med_PPNW);
set(HO2,'Position',[0.755 0.69 0.05 0.04])

leg_line = findobj(HO2,'type','line');
xd2  =get(leg_line(2),'XData');
xd2(1)=xd2(1)+0.025;
xd2(2)=xd2(2)-0.025;
set(leg_line(2),'XData',xd2)
text(735250,kertoja*13.5,'(c) 9 m','Fontsize',font_med_PPNW)

set(gcf,'Color','White')

%========================== CO2 ============================

CO2_mod159 = interp1(zz+dz/2,CO2zt,[1 5 9]);

figure(430)
clf
set(gcf, 'Outerposition',outpos) 
ykonen = subplot(311);
plot(tt,CO2_mod159(1,:)/44.01,'b-','LineWidth',1.5);
hold on
plot(CO2ObsAuto(:,1),CO2ObsAuto(:,3),'k.','LineWidth',1.5); %1.5 m
inx=find(CO2Obs(:,2)==1); 
plot(CO2Obs(inx,1)+datenum(year,1,1),1000/44.01*CO2Obs_korjattu(inx,3),'r+','LineWidth',2);
plot(tt,CO2_eqt/44.01,'k:','LineWidth',1.5);
set(gca,'xlim',[datenum(2013,1,1),datenum(2015,1,1)])
datetick('x','mmm')
set(gca,'xticklabel',[]);

klpois = plot([datenum([2013 1 8]),datenum([2014 12 31])],[320,320],'k','LineWidth',1.5);
klpois2 = plot([datenum([2013 1 8]),datenum([2013 1 8])],[320,340],'k','LineWidth',1.5);
klpois3 = plot([datenum([2014 1 1]),datenum([2014 1 1])],[320,340],'k','LineWidth',1.5);
klpois4 = plot([datenum([2014 12 31]),datenum([2014 12 31])],[320,340],'k','LineWidth',1.5);
set(klpois,'Clipping','off')
set(klpois2,'Clipping','off')
set(klpois3,'Clipping','off')
set(klpois4,'Clipping','off')
hold off
text(datenum([2013 7 1]),345,'Calibration period','HorizontalAlignment','center','FontSize',font_med_PPNW)
text(datenum([2014 7 1]),345,'Validation period','HorizontalAlignment','center','FontSize',font_med_PPNW)

leg=legend('Simul.','Autom. obs.','Man. obs.','\itC_{\rmeq}');
set(leg,'fontsize',font_med_PPNW);
set(leg,'Location','SouthOutside','Orientation','Horizontal');
set(leg,'Box','On');

set(leg,'Position',[0.660 0.865 0.05 0.04])
%Selitteen malliviivan muokkaaminen
leg_line = findobj(leg,'type','line');
xd4 = get(leg_line(4),'XData');
xd4(1) = xd4(1)+0.020;
xd4(2) = xd4(2)-0.020;
set(leg_line(4),'XData',xd4)
xd2  =get(leg_line(2),'XData');
xd2(1)=xd2(1)+0.020;
xd2(2)=xd2(2)-0.020;
set(leg_line(2),'XData',xd2)
xd8  =get(leg_line(8),'XData');
xd8(1)=xd8(1)+0.020;
xd8(2)=xd8(2)-0.020;
set(leg_line(8),'XData',xd8)

text(735250,270,'(a) 1 m','Fontsize',font_med_PPNW)
set(gca,'Ylim',[0 300]);
set(gca,'fontsize',font_med_PPNW);

kakonen = subplot(312);
plot(tt,CO2_mod159(2,:)/44.01,'b-','LineWidth',1.5);
hold on
%plot(kuivameteo(jakso_kuiva(1):jakso_kuiva(2),1),CO2_pinnat_korjattu(:,3),'k-','LineWidth',1.5);
inx=find(CO2Obs(:,2)==5); 
plot(CO2Obs(inx,1)+datenum(year,1,1),1000/44.01*CO2Obs_korjattu(inx,3),'r+','LineWidth',2);
hold off
set(gca,'xlim',[datenum(2013,1,1),datenum(2015,1,1)])
ylabel('CO_{2} (mmol m^{-3})','fontsize',font_med_PPNW)
datetick('x','mmm')
set(gca,'XTickLabel',[])
text(735250,270,'(b) 5 m','Fontsize',font_med_PPNW)
set(gca,'Ylim',[0 300]);
set(gca,'fontsize',font_med_PPNW);

kolomone = subplot(313);
plot(tt,CO2_mod159(3,:)/44.01,'b-','LineWidth',1.5);
hold on
inx=find(CO2Obs(:,2)==9); 
plot(CO2Obs(inx,1)+datenum(year,1,1),1000/44.01*CO2Obs_korjattu(inx,3),'r+','LineWidth',2);
plot(datenum([2012 12 25]),10,'k-','LineWidth',1.5);
plot(datenum([2012 12 25]),10,'k:','LineWidth',1.5);
datetick('x','mmm')
text(735250,400,'(c) 9 m','Fontsize',font_med_PPNW)
set(gca,'Ylim',[0 450]);
set(gca,'fontsize',font_med_PPNW);
set(gca,'xlim',[datenum(2013,1,1),datenum(2015,1,1)])

set(kolomone,'Position',[0.09 0.05 0.88 0.26])
set(kakonen,'Position',[0.09 0.36 0.88 0.26])
set(ykonen,'Position',[0.09 0.67 0.88 0.26])
set(gcf,'Color','w')

%========================== pH ============================

if (~exist('pH_Obs', 'var'))
    [pH_Obs,~]=xlsread('I_O\Kuiva_Obs_CO2.xls','pH');
    pH_dates = pH_Obs(:,1);
end

ph_mod = interp1(zz,10.^(-pHt),[0.5 12]);
ph_mod = -log10(ph_mod);

font_print = 9;
liwi = 1;

figure(260)
clf
set(gcf, 'Outerposition',[500 275 700 400])
plot(tt(1),7,'w.','MarkerSize',8)
hold on
plot(tt,ph_mod(1,:),'k-','LineWidth',liwi);
plot(tt,ph_mod(2,:),'k-.','LineWidth',liwi);
plot(tt(1),7,'w.','MarkerSize',8)
plot(pH_dates,pH_Obs(:,3),'r.','MarkerSize',8);
plot(pH_dates,pH_Obs(:,4),'b.','MarkerSize',8);
plot(pH_dates,pH_Obs(:,2),'.','Color',[0.55 0.27 0.07],'MarkerSize',8);
datetick('x','mmm')
ylabel('pH','fontsize',font_print)
xlim([datenum([2013 1 1]) datenum([2014 1 1])])
ylim([5.5 7.5])
set(gca,'fontsize',font_print)
H11 = legend({'Simulation','surface','bottom','Measurement','surface','bottom','inlet'},'orientation','vertical','location','northeastoutside');
set(H11,'Position',[0.910 0.700 0.05 0.04])

%Selitteen sisällön muokkaaminen:
phleg_text = findobj(H11,'type','text');
% %Selitteen tekstin paikka
set(phleg_text(7),'Position',[0.10 0.9110 0]);
set(phleg_text(6),'Position',[0.28 0.7740 0]);
set(phleg_text(5),'Position',[0.28 0.6370 0]);
set(phleg_text(4),'Position',[0.10 0.5 0]);
set(phleg_text(3),'Position',[0.28 0.3630 0]);
set(phleg_text(2),'Position',[0.28 0.2260 0]);
set(phleg_text(1),'Position',[0.28 0.0890 0]);
%Selitteen malliviivan muokkaaminen
phleg_line = findobj(H11,'type','line');
set(phleg_line(13),'XData',0.90)
set(phleg_line(12),'XData',[0.10 0.22])
set(phleg_line(10),'XData',[0.10 0.22])
set(phleg_line(7),'XData',0.90)
set(phleg_line(5),'XData',0.16)
set(phleg_line(3),'XData',0.16)
set(phleg_line(1),'XData',0.16)

set(H11, 'Box','off')
set(gca,'Position',[0.0700 0.100 0.7400 0.8400])
set(gcf,'Color','w')
