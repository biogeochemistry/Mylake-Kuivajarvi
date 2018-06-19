% Script to run MyLake C for Kuivajärvi
% by PK, last modified 24.04.2018.

global ies80 

if (~exist('Qlambda','var'))
    
    load Qlambda.txt;
    
    %Temperature observations
    depths = [0.2 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 10 12]'; %Havaintosyvyydet [m]
    [numMx,strMx]=xlsread('I_O\Obs_T_Kuivaj_13_14.xls','temp');
    numMx(numMx(:,:)==-999)=NaN;
    [d1,d2]=size(numMx);
    kjtemp=NaN*zeros(d1*length(depths),5);
    hh=0;
    for ii=1:d1
        AAAA=datenum(strMx(ii+1,1),'dd.mm.yyyy');
        [yr,mm,dd]=datevec(AAAA);
        for jj=1:length(depths)
            kjtemp(hh+jj,1:3)=[yr,mm,dd];
            kjtemp(hh+jj,4)=depths(jj);
            kjtemp(hh+jj,5)=numMx(ii,jj);
        end
        hh=hh+length(depths);
    end
    
end

lake='Kuivajärvi'; %Järven nimi
year=2014; %Mallinnuksen aloitusvuosi
m_start= [2013,1,8];%[2014,4,13]; %[2012,5,5]; %Mallinnuksen aloituspäivä [2013,5,11]
m_stop=[2014,12,31];  %Mallinnuksen lopetuspäivä

% year=1982;
% m_start= [1982,1,1];
% m_stop=[1985,12,31];

% year=2082;
% m_start= [2082,1,1];
% m_stop=[2085,12,31];

initfile='I_O\KJ_init_01_13_new.xls'; 
initsheet = 'lake';
parafile='I_O\KJ_para_13_new_joint_50_50.xls';
inputfile='I_O\Input_KJ_2012-2014.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_baseline_av_61_90.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_baseline_av_71_99.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_baseline_av_81_10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_baseline_av_81_10_osahyy.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c2.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c8.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c4.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h2.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h8.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h4.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_c2_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_c4_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_c8_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_h2_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_h4_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_h8_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_m2_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_m4_direct.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_control_81_10_m8_direct.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_c2_direct_discharge.xls';
%%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_c4_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_c8_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_h2_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_h4_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_h8_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m2_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m4_direct_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m8_direct_discharge.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m4_direct_discharge_eivalumamuutosta.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m4_direct_discharge_smoothaamaton.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c4_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c4_discharge_DOC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c4_discharge_DIC.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_c4_direct_discharge_DIC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_c4_direct_discharge_DOC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_h4_direct_discharge_DIC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_h4_direct_discharge_DOC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m4_direct_discharge_DIC.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\input_KJ_scenario_71_00_m4_direct_discharge_DOC.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_c4_direct_discharge_DIC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_c4_direct_discharge_DOC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_h4_direct_discharge_DIC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_h4_direct_discharge_DOC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_m4_direct_discharge_DIC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_m4_direct_discharge_DOC10.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_c4_direct_discharge_DIC20.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_c4_direct_discharge_DOC20.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_h4_direct_discharge_DIC20.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_h4_direct_discharge_DOC20.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_m4_direct_discharge_DIC20.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\direct\DC_changes\input_KJ_scenario_71_00_m4_direct_discharge_DOC20.xls';

%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c2_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c4_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_c8_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h2_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h4_discharge.xls';
%inputfile='C:\MyTemp\DIC\Mylake_v12_package_121107\KJ_application\input_KJ_scenario_71_00_h8_discharge.xls';

tic
[zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,DOCzt1,DOCzt2,DOCzt3,DOCtfrac,...
        Daily_BB1t,Daily_BB2t,Daily_BB3t,Daily_PBt,DICzt,CO2zt,O2zt,O2_sat_relt,O2_sat_abst,POCzt,BODzt,Qzt_sed,lambdazt,...
    P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt,surfaceflux,oxygenflux,CO2_eqt,K0t,O2_eqt,K0_O2t,...
    CO2_ppmt,dO2Chlt,POC1tfrac,dO2SODt,dO2DOCt,pHt,testi3t,T_sedt]...
    = solvemodel_vC(m_start,m_stop,initfile,initsheet,inputfile,'timeseries', parafile,'lake',Qlambda);
run_time=toc

tt_mod = tt - datenum(year,1,1); %time now scaled so that it begins from the 1 january of the "year" (=0)

%=Temperature profile observations (tt_mod, z, T)
TempObs=[datenum(kjtemp(:,1:3)) - datenum(year,1,1),kjtemp(:,4), kjtemp(:,5)];

%=align temperature observations with model results
alku=[1;find(diff(TempObs(:,1))~=0)+1];
loppu=[find(diff(TempObs(:,1))~=0); length(TempObs)];

for i=1:length(alku)
    inxt=find(tt_mod==TempObs(alku(i),1));
    if (isempty(inxt)==0)
        TempMod(alku(i):loppu(i))=interp1(zz,Tzt(:,inxt),TempObs(alku(i):loppu(i),2));
    else
        TempMod(alku(i):loppu(i))=NaN;
    end
end

%==========================================================================

zlim = [0 max(zz)];   %Järven syvyysasteikon rajat
tlim = [min(tt_mod) max(tt_mod)];  %Aika-asteikon rajat

% thermocline depth
zt = MixStat(12,:);

figure(10)   %Jään ja lumen paksuus sekä veden lämpötilaprofiili
clf
subplot(4,1,2:4) %Veden lämpötilaprofiili
pcolor(tt_mod,zz,Tzt)
shading interp
axis ij
hold on
plot(tt_mod,zt,'k','LineWidth',1);   %Termokliinin syvyyskäyrä
hold off
datetick('x','mmm'); %Vaaka-akselin jakoviivojen päivämääräotsikon muoto
set(gca,'ylim',zlim); %Pystyakselin rajat
caxis([0 25]);   %Lämpötila-asteikon rajat
colorbar;
set(gca,'fontsize',9);
ylabel('Depth (m)')   %Pystyakselin otsikko
set(gca,'TickDir','out')

subplot(4,1,1);   %Jään ja lumen paksuus
plot(tt_mod,His(1,:)+His(2,:),'-r',tt_mod,His(1,:)-His(3,:),':c',tt_mod,His(1,:),'-b')
% Jään ja lumen yhteispaksuus (punainen viiva), teräsjään paksuus
% (sinivihreä pisteviiva), kokonaisjäänpaksuus (sininen viiva)

%H(2)=H(2)+0.23; %0.226
%set(gca,'Position',H);
set(gca,'fontsize',9);
set(gca,'ylim',[0 1]);   %Pystyakselin rajat
datetick('x','mmm');   %Vaaka-akselin jakoviivojen päivämääräotsikon muoto
set(gca,'XTickLabel',[]);
set(gca,'YTick',[0.2 0.4 0.6 0.8 1]); %Pystyakselin jaotus
set(gca,'TickDir','out')
ylabel('Hice, Hsnow (m)');   %Pystyakselin otsikko

grid on;

dz=zz(2)-zz(1);
depths_f18 = [1 5 10];
Tzt_f18=interp1(zz+dz/2,Tzt,depths_f18);
tlims=[datenum(m_start) datenum(m_stop)];

figure(18) %Mallinnetut ja havaitut päivittäiset lämpötilat eri kerroksissa
clf
subplot(311) %Pintakerros, 0,5 metriä - T-havainnot vuosina 90-99 ja 04-07 syvyydeltä 0 m
inx=find(TempObs(:,2)==depths_f18(1));
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2); %Havainnot
hold on
plot(tt_mod+datenum(year,1,1),Tzt_f18(1,:),'-'); %Mallinnettu lämpötila 0,5 metrissä
set(gca,'ylim',[0 25]); %Lämpötila-asteikon rajat
ylabel('^oC','fontsize',9) %Pystyakselin otsikko
title(['Temperature ',num2str(depths_f18(1)),' m'],'fontweight','bold') %Kuvaajan otsikko
datetick('x','mm yy') %Vaaka-akselin jakoviivojen päivämääräotsikon muoto
grid on
set(gca,'fontsize',9)
set(gca,'xticklabel',[]);
set(gca,'xlim',tlims) %Aika-akselin rajat

subplot(312) %Lämpötilan harppauskerros, 3 metriä (b=c-1)
inx=find(TempObs(:,2)==depths_f18(2)); %havainnot 2000-2003
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt_f18(2,:),'-'); %Mallinnettu lämpötila 2 metrissä
set(gca,'ylim',[0 25]); %Lämpötila-asteikon rajat
ylabel('^oC','fontsize',9) %Pystyakselin otsikko
title(['Temperature ',num2str(depths_f18(2)),' m'],'fontweight','bold') %Kuvaajan otsikko
datetick('x','mm yy') %Vaaka-akselin jakoviivojen päivämääräotsikon muoto
grid on
set(gca,'fontsize',9)
set(gca,'xticklabel',[]);
set(gca,'xlim',tlims) %Aika-akselin rajat

subplot(313) %Alusvesi, 5 metriä - T-havainnot vuosina 90-99 ja 04-07 syvyydeltä 5 m
inx=find(TempObs(:,2)==depths_f18(3)); %havainnot 2000-2003
plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r+','linewidth',2);
hold on
plot(tt_mod+datenum(year,1,1),Tzt_f18(3,:),'-');
set(gca,'ylim',[0 25]); %Lämpötila-asteikon rajat
ylabel('^oC','fontsize',9) %Pystyakselin otsikko
xlabel('year','fontsize',9) %Vaaka-akselin otsikko
title(['Temperature ',num2str(depths_f18(3)),' m'],'fontweight','bold') %Kuvaajan otsikko
datetick('x','mm yy') %Vaaka-akselin jakoviivojen päivämääräotsikon muoto
grid on
set(gca,'fontsize',9)
set(gca,'xlim',tlims) %Aika-akselin rajat

%===============
man_figures; 
evaluation_KJ_2013;

ycp = 1.2698;
C_Chla = 50/ycp;
F_ALPOC = 1000;
f_POC_1 = POC1tfrac;
f_POC_2 = squeeze(testi3t(:,:,6));
P_POCzt = f_POC_1.*POCzt/(C_Chla*ycp)+f_POC_2.*POCzt/F_ALPOC;
Ptot = Pzt+PPzt+DOPzt+(Chlzt+Czt)./ycp+P_POCzt;
Chltot = Chlzt+Czt;
P_Chlzt = Chltot./ycp;
