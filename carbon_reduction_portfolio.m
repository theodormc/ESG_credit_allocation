% carbon neutral portfolio
% first compute the carbon intensity of each company 
% depending on the portfolio structure
addpath('C:\Users\XYZW\Documents\Matlab dissertation - new version\Functions')
addpath('C:\Users\XYZW\Documents\Matlab dissertation - new version\Markowitz_problems\risk_contributions_budget')
CE1 = [75,5000,720,50,2500,25,30000,5];
CE2 = [75,5000,1030,350,4500,5,2000,64];
CE3 = [24000,15000,1210,550,500,187,30000,199];
Y = [300,320,125,100,200,102,100,25];
b = [0.23,0.19,0.17,0.13,0.09,0.08,0.06,0.05];
W = 1;MktCap_total = 2500;
carb_intens1 = CE1./Y;carb_intens2 = CE2./Y;carb_intens3 = CE3./Y;
carb_intens12 = (CE1+CE2)./Y;carb_intens123 = (CE1+CE2+CE3)./Y;
mkt_caps = b*MktCap_total;
CE1_func = @(x) W*x./mkt_caps.*CE1;
CE2_func = @(x) W*x./mkt_caps.*CE2;
CE3_func = @(x) W*x./mkt_caps.*CE3;
Y_func = @(x) W*x./mkt_caps.*Y;
CI1_portf = CE1_func(b)/Y_func(b);CI2_portf = CE2_func(b)/Y_func(b);
CI3_portf = CE3_func(b)/Y_func(b);
%%
sigs = [0.22,0.2,0.25,0.18,0.35,0.23,0.13,0.29];
carb_intens_12 = carb_intens1+carb_intens2;
carb_intens_123 = carb_intens1+carb_intens2+carb_intens3;
ro = [1,0.8,0.7,0.6,0.7,0.5,0.7,0.6;0.8,1,0.75,0.65,0.5,0.6,0.5,0.65;...
    0.7,0.75,1,0.8,0.7,0.7,0.7,0.7;0.6,0.65,0.8,1,0.85,0.8,0.75,0.75;...
    0.7,0.5,0.7,0.85,1,0.6,0.8,0.65;0.5,0.6,0.7,0.8,0.6,1,0.5,0.7;0.7,0.5,0.7,0.75,0.8,0.5,1,0.8;...
    0.6,0.65,0.7,0.75,0.65,0.7,0.8,1];
CI12_portf = (CE1_func(b)+CE2_func(b))/Y_func(b);
cov_mat = (sigs'.*sigs).*ro;
func = @(x) (x-b)*cov_mat*(x-b)';
n = length(sigs);R = 0.2;
C12 = carb_intens_12;D12 = (1-R)*b*carb_intens_12';
C123 = carb_intens_123;D123 = (1-R)*b*carb_intens_123';
C_HCIS_12 = [C12;0,-1,0,0,-1,0,-1,0];
D_HCIS_12 = [D12;-0.34];
C_HCIS_123 = [C123;0,-1,0,0,-1,0,-1,0];
D_HCIS_123 = [D123;-0.34];
%%
str12 = fmincon(func,ones(1,n)/n,C12,D12,ones(1,n),1,zeros(1,n),ones(1,n));
str123 = fmincon(func,ones(1,n)/n,C123,D123,ones(1,n),1,zeros(1,n),ones(1,n));
str12_HCIS = fmincon(func,ones(1,n)/n,C_HCIS_12,D_HCIS_12,ones(1,n),1,zeros(1,n),ones(1,n));
str123_HCIS = fmincon(func,ones(1,n)/n,C_HCIS_123,D_HCIS_123,ones(1,n),1,zeros(1,n),ones(1,n));
Carb_em_12 = (CE1+CE2)./mkt_caps;D_em_12 = b*(1-R)*((CE1+CE2)./mkt_caps)';
str12_em = fmincon(func,ones(1,n)/n,Carb_em_12,D_em_12,ones(1,n),1,zeros(1,n),ones(1,n));
Carb_em_12_HCIS = [Carb_em_12;0,-1,0,0,-1,0,-1,0];
D_em_12_HCIS = [D_em_12;-0.34];
str12_em_HICS = fmincon(func,ones(1,n)/n,Carb_em_12_HCIS,D_em_12_HCIS,ones(1,n),1,zeros(1,n),ones(1,n));
tracking_err_vol = func(str12);
carb_intensities_abs = [100.5,57.2,250.4,352.3,27.1,54.2,78.6,426.7];
str_example = min_risk_decarb2(sigs,ro,b,carb_intensities_abs,0.2);
str_example2 = min_risk_decarb2(sigs,ro,b,carb_intensities_abs,0.3);
%%
% order-statistic approach
% eliminate the worst carbonized portfolio
pos = find(CE1+CE2+CE3==max(CE1+CE2+CE3));
z = zeros(1,n);z(pos)=1;
Aeq = [ones(1,n);z];beq = [1;0];
str12_wo_worst = fmincon(func,ones(1,n)/n,C12,D12,Aeq,beq,zeros(1,n),ones(1,n));
sorted_carb_intens = sort(CE1+CE2+CE3);
pos_123 = arrayfun(@(i) find(CE1+CE2+CE3==sorted_carb_intens(i)),1:3);
z=zeros(3,n);z(1,pos_123(1))=1;z(2,pos_123(2))=1;z(3,pos_123(3))=1;
z2 = zeros(2,n);z2(1,pos_123(1))=1;z2(2,pos_123(2))=1;
Aeq = [ones(1,n);z];beq = [1;zeros(3,1)];
Aeq2 = [ones(1,n);z2];beq2 = [1;zeros(2,1)];
str12_wo_worst3 = fmincon(func,ones(1,n)/n,C12,D12,Aeq,beq,zeros(1,n),ones(1,n));
str123_wo_worst3 = fmincon(func,ones(1,n)/n,C123,D123,Aeq,beq,zeros(1,n),ones(1,n));
str12_wo_worst2 = fmincon(func,ones(1,n)/n,C12,D12,Aeq2,beq2,zeros(1,n),ones(1,n));
str123_wo_worst2 = fmincon(func,ones(1,n)/n,C123,D123,Aeq2,beq2,zeros(1,n),ones(1,n));
%%
clc
str_alloc_order_stat = order_stat_CA(CE1,CE2,CE3,b,Y,R,sigs,ro);
str_func_CI_decarb_R12 = min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,[]);
str_func_CI_decarb_R123 = min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,[]);
str_func_CI_decarb_R12_abs = min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,[],'type','absolute');
str_func_CI_decarb_R123_abs = min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,[],'type','absolute');
str_func_CI_decarb_R12_emission = min_risk_decarb(sigs,ro,b,CE1+CE2,Y,0.2,[],'type','relative','carbon','WACI');
str_func_CI_decarb_R123_emission = min_risk_decarb(sigs,ro,b,CE1+CE2+CE3,Y,0.2,[],'type','relative','carbon','WACI');
%%
function str = min_risk_decarb2(sigs,corr,bench,CIs,R)
   cov_mat = sigs'.*sigs.*corr;
   func = @(x) (x-bench)*cov_mat*(x-bench)';
   n = length(sigs);
   C = CIs;D = (1-R)*bench*CIs';
   str = fmincon(func,ones(1,n)/n,C,D,ones(1,n),1,zeros(1,n),ones(1,n));
end
%%
