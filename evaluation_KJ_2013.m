%Performance metrics for T, O2, CO2 in 2013

%NOTE! DO in mmol/l

% Tästä observations used in calibration =======================================

%CO2 ------------------
%if (exist('hiilidhavainnot', 'var') == 0)
    [hiilidhavainnot,StrMx] = xlsread...
        ('I_O\Kuiva_Obs_CO2.xls','CO2MCMC');
    dumda2 = datevec(StrMx(:,1),'dd.mm.yyyy');
    hiilidhavainnot(1,:) = [];
    %Poistetaan kokonaan vuoden -14 havainnot
    dumda2(47:end,:) = [];
    hiilidhavainnot(47:end,:) = [];
    %Perättäisten päivien havainnot pois
    dumda2(11:14,:) = [];
    hiilidhavainnot(11:14,:) = [];
    %Myös 7.1.13 pois, koska sille ei ole simulaatiota!
    dumda2(1,:) = [];
    hiilidhavainnot(1,:) = [];
%end

KJ_obs_CO2 = [dumda2(:,1:3) hiilidhavainnot];

%DO ---------------------------

%if (~exist('happihavainnot_eval', 'var'))
    [happihavainnot_eval,~] = xlsread...
        ('I_O\KJ_Obs_O2.xls','oxygen_MCMC');
    happihavainnot_eval(1,:) = [];
    paivavek = happihavainnot_eval(:,1:3);
    
    %Vuosi -14 pois
    happihavainnot_eval(47:end,:) = [];
    paivavek(47:end,:) = [];
    %Perättäisten päivien havainnot pois
    happihavainnot_eval(11:14,:) = [];
    paivavek(11:14,:) = [];
    %Myös 7.1.13 pois, koska sille ei ole simulaatiota!
    paivavek(1,:) = [];
    happihavainnot_eval(1,:) = [];
%end

KJ_obs_O2 = happihavainnot_eval;

if(any(datenum(dumda2)-datenum(paivavek))~=0)
    keyboard
    disp('Happi- ja hiilidioksidihavaintopäivät eivät ole samat!')
end
%---------------------

Tobs=[happihavainnot_eval hiilidhavainnot]; %happi ja hiilidioksidi samassa

happihavainnot_eval_13 = Tobs(:,1:26);
CO2havainnot_13 = Tobs(:,27:end);
CO2havaintopaivat_13 = Tobs(:,1:3);

%Open water season manual CO2 observations
%CO2havainnot_13 = Tobs(9:39,27:end);
%CO2havaintopaivat_13 = Tobs(9:39,1:3);

%Automatic manual CO2 observations during simulated open water season
%2.5.-24.11.2013
CO2_auto_havainnot_13 = 44.01/1000*CO2ObsAuto(115:321,3);
CO2_auto_havaintopaivat_13 = (735356:735562)';

%T ==========================================================================

thin_switch = 0;
[numMx,strMx]=xlsread('I_O\Obs_T_Kuivaj_13_14.xls','temp');
[d1,d2]=size(numMx);
Obs_T_02=NaN*zeros(d1*length(depths),5);
if(thin_switch==1)    
   hh=0;
   for ii=1:5:d1
      AAAA=datenum(strMx(ii+1,1),'dd.mm.yyyy');
      [yr,mm,dd]=datevec(AAAA);
      for jj=1:length(depths)
         Obs_T_02(hh+jj,1:3)=[yr,mm,dd];
         Obs_T_02(hh+jj,4)=depths(jj);
         Obs_T_02(hh+jj,5)=numMx(ii,jj);
      end
      hh=hh+length(depths);
   end
else
   hh=0;
   for ii=1:d1
      AAAA=datenum(strMx(ii+1,1),'dd.mm.yyyy');
      [yr,mm,dd]=datevec(AAAA);
      for jj=1:length(depths)
         Obs_T_02(hh+jj,1:3)=[yr,mm,dd];
         Obs_T_02(hh+jj,4)=depths(jj);
         Obs_T_02(hh+jj,5)=numMx(ii,jj);
      end
      hh=hh+length(depths);
   end
end

nans = find(isnan(Obs_T_02(:,1)));
Obs_T_02(nans,:) = [];

[yr,mm,dd]=datevec(datenum(strMx(2:end,1),'dd.mm.yyyy'));
Obs_T_02_calval=[yr,mm,dd,numMx(:,:)];

Obs_T_02_calval(1:7,:) = [];

nans = find(isnan(Obs_T_02_calval(:,5)));
Obs_T_02_calval(nans,:) = [];
%Year 2013
Obs_T_02_calval(357:end,:) = [];

%Observed open water season
%Obs_T_02_calval(322:end,:) = [];
%Obs_T_02_calval(1:113,:) = [];
%Simulated open water season
%Obs_T_02_calval(320:end,:) = [];
%Obs_T_02_calval(1:114,:) = [];

Tobs_calval=Obs_T_02_calval;

lukumaara = 1;
iii = lukumaara;

T_summat = NaN*zeros(lukumaara,16);
r2s = NaN*zeros(lukumaara, 16);
rmses = NaN*zeros(lukumaara, 16);
NSEs = NaN*zeros(lukumaara, 16);
Pbiass = NaN*zeros(lukumaara, 16);
Nbiass = NaN*zeros(lukumaara, 16);
rmsenus = NaN*zeros(lukumaara, 16);
Biass = NaN*zeros(lukumaara, 16);

O2_summat = NaN*zeros(lukumaara,23);
r2sO2 = NaN*zeros(lukumaara, 23);
rmsesO2 = NaN*zeros(lukumaara, 23);
NSEsO2 = NaN*zeros(lukumaara, 23);
PbiassO2 = NaN*zeros(lukumaara, 23);
NbiassO2 = NaN*zeros(lukumaara, 23);
rmsenusO2 = NaN*zeros(lukumaara, 23);
BiassO2 = NaN*zeros(lukumaara, 23);

CO2_summat = NaN*zeros(lukumaara,8);
r2sCO2 = NaN*zeros(lukumaara, 8);
rmsesCO2 = NaN*zeros(lukumaara, 8);
NSEsCO2 = NaN*zeros(lukumaara, 8);
PbiassCO2 = NaN*zeros(lukumaara, 8);
NbiassCO2 = NaN*zeros(lukumaara, 8);
rmsenusCO2 = NaN*zeros(lukumaara, 8);
BiassCO2 = NaN*zeros(lukumaara, 8);

CO2auto_summat = NaN*zeros(lukumaara,1);
r2sCO2auto = NaN*zeros(lukumaara, 1);
rmsesCO2auto = NaN*zeros(lukumaara, 1);
NSEsCO2auto = NaN*zeros(lukumaara, 1);
PbiassCO2auto = NaN*zeros(lukumaara, 1);
NbiassCO2auto = NaN*zeros(lukumaara, 1);
rmsenusCO2auto = NaN*zeros(lukumaara, 1);
BiassCO2auto = NaN*zeros(lukumaara, 1);

tt_mod = tt - datenum(year,1,1);
Obs_T_tt_mod=datenum(Tobs_calval(:,1:3))- datenum(year,1,1);

%Tfit_all=interp1(zz+dz/2,Tzt,[1;5;10])'; %interpolated to observed depth level, dim(t, z); requires dz<2*z_obs(1)
Tfit_all=interp1(zz+dz/2,Tzt,[0.25;0.5;1;1.5;2;2.5;3;3.5;4;4.5;5;6;7;8;10;12])';
Tfit_sub=interp1(tt_mod,Tfit_all,Obs_T_tt_mod); %dim(t, z)

Out_zlevel = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
SS_pure_T_new=nansum((Tfit_sub - Tobs_calval(:, 3+Out_zlevel)).^2) %Sum of squares
Nobs_T=sum(isnan(Tfit_sub - Tobs_calval(:, 3+Out_zlevel))==0) %Number of observations in the model period

T_summat(iii,:) = SS_pure_T_new;

[r2s(iii,:),rmses(iii,:),NSEs(iii,:),Pbiass(iii,:),Nbiass(iii,:),rmsenus(iii,:),Biass(iii,:)] = evalstat(Tobs_calval(:, 3+Out_zlevel),Tfit_sub);
 
Obs_O2_tt_mod=datenum(happihavainnot_eval_13(:,1:3))- datenum(year,1,1); %modified timevector for observations
%0 m -> 0.25 m
O2fit_all=interp1(zz+dz/2,O2zt,[0.25 0.5:0.5:9 10 11 12 12.5]')'; %O2 interpolated to observed depth level, dim(t, z); requires dz<2*z_obs(1)
O2fit_sub=interp1(tt_mod,O2fit_all,Obs_O2_tt_mod); %dim(t, z)

Obs_CO2_tt_mod=datenum(CO2havaintopaivat_13(:,1:3)) - datenum(year,1,1);
CO2fit_all=interp1(zz+dz/2,CO2zt,[0.25 1 3 5 7 9 11 12]')';
CO2fit_sub=interp1(tt_mod,CO2fit_all,Obs_CO2_tt_mod);

SS_pure_O2_new=nansum((0.001*O2fit_sub(:,1:23) - happihavainnot_eval_13(:, 4:26)).^2) %Sum of squares
SS_pure_CO2_new=nansum((0.001*CO2fit_sub(:,1:8) - CO2havainnot_13(:, 1:8)).^2)
Nobs_O2=sum(isnan(0.001*O2fit_sub(:,1:23) - happihavainnot_eval_13(:, 4:26))==0)
Nobs_CO2=sum(isnan(0.001*CO2fit_sub(:,1:8) - CO2havainnot_13(:, 1:8))==0)

O2_summat(iii,:) = SS_pure_O2_new;
CO2_summat(iii,:) = SS_pure_CO2_new;

[r2sO2(iii,:),rmsesO2(iii,:),NSEsO2(iii,:),PbiassO2(iii,:),NbiassO2(iii,:),rmsenusO2(iii,:),BiassO2(iii,:)]...
     = evalstat(1000*1/32*happihavainnot_eval_13(:, 4:26),1/32*O2fit_sub(:,1:23));

[r2sCO2(iii,:),rmsesCO2(iii,:),NSEsCO2(iii,:),PbiassCO2(iii,:),NbiassCO2(iii,:),rmsenusCO2(iii,:),BiassCO2(iii,:)]...
     = evalstat(CO2havainnot_13(:, 1:8),0.001*CO2fit_sub(:,1:8));

eval_T_cal13 = [T_summat; r2s; rmses; NSEs; Pbiass; Nbiass; rmsenus; Nobs_T]; 
eval_O2_cal13 = [O2_summat; r2sO2; rmsesO2; NSEsO2; PbiassO2; NbiassO2; rmsenusO2; Nobs_O2];
eval_CO2_cal13 = [CO2_summat; r2sCO2; rmsesCO2; NSEsCO2; PbiassCO2; NbiassCO2; rmsenusCO2; Nobs_CO2];


Obs_CO2auto_tt_mod=CO2_auto_havaintopaivat_13 - datenum(year,1,1);
CO2autofit_all=interp1(zz+dz/2,CO2zt,1.5)';
CO2autofit_sub=interp1(tt_mod,CO2autofit_all,Obs_CO2auto_tt_mod);

SS_pure_CO2auto_new=nansum((0.001*CO2autofit_sub(:,1) - CO2_auto_havainnot_13(:, 1)).^2)

Nobs_CO2auto=sum(isnan(0.001*CO2autofit_sub(:,1) - CO2_auto_havainnot_13(:, 1))==0)

CO2auto_summat(iii,:) = SS_pure_CO2auto_new;

[r2sCO2auto(iii,:),rmsesCO2auto(iii,:),NSEsCO2auto(iii,:),PbiassCO2auto(iii,:),NbiassCO2auto(iii,:),rmsenusCO2auto(iii,:),BiassCO2auto(iii,:)]...
     = evalstat(CO2_auto_havainnot_13(:, 1),0.001*CO2autofit_sub(:,1));

eval_CO2auto_cal13 = [CO2auto_summat; r2sCO2auto; rmsesCO2auto; NSEsCO2auto; PbiassCO2auto; NbiassCO2auto; rmsenusCO2auto; Nobs_CO2auto];

%{
tg_T_cal = [rmsenus;Nbiass];
tg_O2_cal = [rmsenusO2;NbiassO2];
tg_CO2_cal = [rmsenusCO2;NbiassCO2];

save tgcalKJ.mat tg_T_cal  tg_O2_cal tg_CO2_cal
save evalcalKJ.mat eval_T_cal13 eval_O2_cal13 eval_CO2_cal13
%}

rs_T = zeros(16,2,2);
ps_T = zeros(16,2,2);
for g = 1:16
    [rs_T(g,:,:),ps_T(g,:,:)] = corrcoef(Tobs_calval(:, 3+g),Tfit_sub(:,g));
end

final_R2s_T = squeeze(rs_T(:,2,1)).^2;
final_ps_T = squeeze(ps_T(:,2,1));

rs_O = zeros(23,2,2);
ps_O = zeros(23,2,2);
for g = 1:23
    einanit = find(~isnan(happihavainnot_eval_13(:,3+g)));
    [rs_O(g,:,:),ps_O(g,:,:)] = corrcoef(happihavainnot_eval_13(einanit, 3+g),0.001*O2fit_sub(einanit,g));
end
final_R2s_O = squeeze(rs_O(:,2,1)).^2;
final_ps_O = squeeze(ps_O(:,2,1));

rs_C = zeros(8,2,2);
ps_C = zeros(8,2,2);
for g = 1:8
    einanit = find(~isnan(CO2havainnot_13(:,g)));
    [rs_C(g,:,:),ps_C(g,:,:)] = corrcoef(CO2havainnot_13(einanit, g),0.001*CO2fit_sub(einanit,g));
end
final_R2s_C = squeeze(rs_C(:,2,1)).^2;
final_ps_C = squeeze(ps_C(:,2,1));
