%利用压缩感知实现南京明城墙仪凤门地区层析成像(共26景，TerraSAR），预处理使用了GAMMA软件，具体过程参考GAMMA_script.txt.
%截取了原图像的一块45*45的仪凤门区域。
%影像集In以20150206（第4个文件夹，距离基线为中间）为主影像，其它影像经过了配准（offset_pwr）。
%尤江彬
%中国科学院遥感与数字地球研究所，radi
%start：2017.2.8

clc;clear;close all;
%% 参数设置&初始化
% c = 3e8;    %光速
wave_length = 0.031;  %波长
incidence_angle = 37.3*pi/180;   %入射角
j = sqrt(-1);   %虚数单位
center_slant_range = 633773.4808;   %中心处斜距
near_slant_range = 633753.4739;   %近距处斜距
far_slant_range = 633793.4877;   %远距处斜距
baseline_number = 26;
Tab_used = ones(1,baseline_number);
img_width = 45; %range,column
img_height = 45; %azimuth,row
% bn=[-296.3659,-229.8177,-228.5793,-200.2995,-176.9877,-161.6247,-126.1341,-107.9265,-104.9532,-86.2247,-73.6160,-71.7670,-56.6469,-51.1142,-38.2620,-9.2436,0,16.5665,36.4709,51.5930,68.1936,74.6618,117.4378,177.7292,227.8118,300.0834];
% bp=[38.3175, 302.8135,87.4850,384.1936,386.7169,343.0427,144.1888,311.3235,131.6339,185.6147,417.5037,505.5218,129.4644,405.0742,457.4194,20.0097,0,278.5639,51.3897,101.3265,552.2287,183.2972,83.8275,133.8108,220.9207,152.3898];
%上面是针对整幅影像计算的基线；下面是根据仪凤门段影像计算得到的基线
bn = [-294.9578,-230.1109,-227.7124,-201.7206,-178.5367,-162.7868,-125.6034,-108.7292,-104.2799,-85.9184,-75.1820,-73.8293,-56.1035,-52.7090,-40.6110,-7.8088,0,16.0799,37.7794,52.3234,65.2059,74.4328,119.3556,178.1374,228.0251,300.7758];
bp = [34.4957,299.1064,84.0729,381.0866,383.6235,340.1932,141.3058,308.9355,128.9492,183.0384,415.1938,503.4068,126.8287,402.8447,455.1216,17.2105,0,276.6002,48.7611,98.9504,550.5676,181.5631,81.1736,131.9808,219.6773,151.1835];
bn = bn.*Tab_used;
bp = bp.*Tab_used;
% figure;
% bar(bn);   %基线分布概览
date = [20140723,20141202,20150115,20131123,20130827,20140106,20140518,20131101,20140814,20140701,20140128,20131215,20140927,20130918,20131010,20130805,20150206,20141110,20130531,20140905,20140219,20141019,20141224,20130714,20130622,20140609];
%借助index实现由时间基线顺序向空间基线顺序的转换
bn_t=[-26.0646,163.5846,114.1326,-71.5183,-243.2902,-117.5589,-105.4921,-173.3177,-266.5042,-139.0797,-227.4289,-140.1063,0,-189.7014,236.5578,-150.1771,-358.7338,-168.3412,-11.6086,-120.1253,10.3680,-48.4379,-294.8236,55.1146,-291.6025,-63.1099];
bn_t_20140128_cfl=[113.83,303.16,254.63,68.81,-102.96,23.02,35.22,-33.09,-125.49,0.66,-87.22,0,140.29,-49.15,375.85,-10.11,-219.30,-28.35,127.44,19.78,149.93,91.43,-154.36,195.37,-151.47,76.89];%陈老师的垂直基线数据，经检验与上面我自己的是一致的
[sbn,index] = sort(bn_t); 
b_t_t = [-1.688,-1.627,-1.567,-1.507,-1.447,-1.386,-1.326,-1.266,-1.205,-1.145,-1.085,-1.025,-0.964,-0.723,-0.663,-0.603,-0.542,-0.482,-0.422,-0.362,-0.301,-0.241,-0.181,-0.121,-0.060,0.000];
%b_k_t1 = [18.9,27.3,30.1,32,27.2,24.8,22.8,17.1,14.4,8.5,7.1,7.7,2.1,20.3,25.4,22.6,31.2,22.4,24.5,24.8,21,13.1,2.2,6,6.3,3.2];
b_k_t = [24,31,32,36,32,29,28,20,18,12,13,13,5,25,31,25,35,24,29,28,27,17,6,12,10,9];%陈老师的温度数据，与我的还是差别挺大的。。。
b_h_t = [14.464,16.799,19.88,21.335,19.317,16.212,15.696,11.975,9.301,5.653,5.242,5.369,4.603,14.291,17.24,18.54,24.337,19.126,18.635,18.724,13.966,8.243,3.221,3.932,5.644,3.446];
b_t = b_t_t(index);
b_k = b_k_t(index);
b_k = b_k - b_k(17);
b_h = b_h_t(index);
b_h = b_h - b_h(17);

bap = max(bn)-min(bn);  %cross_range孔径大小（基线孔径）
ro = 0.5*wave_length*center_slant_range/bap;
dv = 1;  %cross_range向分辨率
v = -20:dv:100;       %cross_range采样点
r = linspace(near_slant_range,far_slant_range,img_width);   %斜距分布
range_spacing = (far_slant_range-near_slant_range)/img_width;      %影像斜距向分辨率
azimuth_spacing = 2.04;  %影像方位向分辨率
%% 
[base_filename,base_filepath] = uigetfile('*.*','请选择基线文件');% 基线文件应包括两列数据：1.日期 2.垂直基线
str_base_list = [base_filepath,base_filename];
baseline_file=load(str_base_list);
master=baseline_file(1,2);%master应由系统前面的功能确定
file_name_date=baseline_file(:,3);
file_name_date=sort(file_name_date);
baseline_number=length(file_name_date);
master_par=strcat(path,'piece/',num2str(master),'.rslc.par');%master的par文件同样应由前面搞定
header=fopen(master_par);
line_number=0;
C=299792458;
while 1     %对主影像的par文件进行参数抓取
    line_number=line_number+1;
    line_content=fgetl(header);
    if line_number==11
        line_content(find(isspace(line_content)))=[];
        img_width=str2num(cell2mat(regexp(line_content,'\d','match')));
    end
    if line_number==12
        line_content(find(isspace(line_content)))=[];
        img_height=str2num(cell2mat(regexp(line_content,'\d','match')));
    end
    if line_number==22
        line_content(find(isspace(line_content)))=[];
        range_spacing=str2num(line_content(length('range_pixel_spacing:')+1:length(line_content)-1));
    end
    if line_number==23
        line_content(find(isspace(line_content)))=[];
        azimuth_spacing=str2num(line_content(length('azimuth_pixel_spacing:')+1:length(line_content)-1));
    end
    if line_number==24
        line_content(find(isspace(line_content)))=[];
        near_slant_range=str2num(line_content(length('near_range_slc:')+1:length(line_content)-1));
    end
    if line_number==30
        line_content(find(isspace(line_content)))=[];
        incidence_angle=str2num(line_content(length('incidence_angle:')+1:length(line_content)-7));
        incidence_angle = incidence_angle*pi/180;
    end
    if line_number==33
        line_content(find(isspace(line_content)))=[];
        wave_length=C/str2num(line_content(length('radar_frequency:')+1:length(line_content)-2));
        break;
    end
end
%温度基线
temperature=zeros(baseline_number,1);
temperature_data=importdata('/home/mapeifeng/mount/E/former_Ubantu/100/HKO_noheader.txt');
for index=1:baseline_number
    day=temperature_data(find(temperature_data(:,1)==file_name_date(index)),:);
    %temperature(index)=day(find(day(:,2)==10),7);
    %temperature(index)=(temperature(index)+day(find(day(:,2)==11),7))/2;
    temperature(index)=day(1,3); %暂定温度数据文件仅由三列组成：1.日期 2时刻 3.温度
    if file_name_date(index)==master
        master_temperature=temperature(index);
    end
end
temperature(:)=temperature(:)-master_temperature;%slave-master
%空间基线
% time=baseline_file(:,5);
% sort_time=sort(time);
% total_mean_baseline=zeros(baseline_number,1);
% for row=1:baseline_number
%     total_mean_baseline(row,1)=baseline_file(find(time==sort_time(row)),4);%read baseline file
% end
% total_mean_baseline=total_mean_baseline';
mean_baseline=zeros(baseline_number,1);
for row=1:baseline_number
    mean_baseline(row,1)=baseline_file(find(baseline_file(:,3)==file_name_date(row)),4);
end
%时间基线
time_baseline=zeros(baseline_number,1);
for row=1:baseline_number
    time_baseline(row,1)=(datenum(num2str(file_name_date(row)),'yyyymmdd')-datenum(num2str(master),'yyyymmdd'))/365;
end
%% SLC文件读取
In = single(zeros(img_height,img_width,baseline_number));   %SLC数据
In_r = single(zeros(img_height,img_width,baseline_number));
In_i = single(zeros(img_height,img_width,baseline_number));
for i = 1:baseline_number
    fid = fopen(['E:\data\yifengmen\new_reg\section',num2str(i),'.slc']);
    In_r(:,:,i) = fread(fid,[img_height,img_width],'*short',2,'ieee-be').';
    frewind(fid);
    fread(fid,1,'*short');
    In_i(:,:,i) = fread(fid,[img_height,img_width],'*short',2,'ieee-be').';
    fclose(fid);
end
fid = fopen('E:\data\TomoSAR\yifengmen\In_r.bin','w');
fwrite(fid,In_r,'int16');
fclose(fid);
fid = fopen('E:\data\TomoSAR\yifengmen\In_i.bin','w');
fwrite(fid,In_i,'int16');
fclose(fid);
%读取
disp('Reading SLCs...');
fid1 = fopen('E:\data\TomoSAR\yifengmen\In_r.bin');
fid2 = fopen('E:\data\TomoSAR\yifengmen\In_i.bin');
for i = 1:baseline_number
    In(:,:,i) = fread(fid1,[img_height,img_width],'int16=>single');
    In(:,:,i) =In(:,:,i) + fread(fid2,[img_height,img_width],'int16=>single')*j;
end
fclose(fid2);
fclose(fid1);

%% 去地平(相位存于文件ph_flt.bin中)
disp('Processing phase flatten...');
ph_flt = single(zeros(img_height,img_width,baseline_number));
for i = 1:baseline_number
    if i == 17||Tab_used(i) == 0
        continue;
    else
        fid = fopen(['E:\data\yifengmen\new_reg\interf\17_',num2str(i),'.int']);
        igram_r = fread(fid,[img_height,img_width],'*single',4,'ieee-be').';
        frewind(fid);
        fread(fid,1,'*single');
        igram_i = fread(fid,[img_height,img_width],'*single',4,'ieee-be').';
        fclose(fid);
        fid = fopen(['E:\data\yifengmen\new_reg\interf\17_',num2str(i),'.flt']);
        iflt_r = fread(fid,[img_height,img_width],'*single',4,'ieee-be').';
        frewind(fid);
        fread(fid,1,'*single');
        iflt_i = fread(fid,[img_height,img_width],'*single',4,'ieee-be').';
        fclose(fid);
        ph_flt(:,:,i) = angle((igram_r+igram_i*j).*(iflt_r-iflt_i*j));
    end
end
In = In.*exp(j*ph_flt);
%% 振幅离差指数阈值法筛选PS点
abs_In = abs(In);
abs_In_2 = abs_In.^2;
avg_In = mean(abs_In,3);
avg_In_2 = mean(abs_In_2,3);
adi = sqrt(avg_In_2-avg_In.^2)./avg_In;
SPSC_index = find(avg_In>100);
PPSC_index = find(adi<0.23);
%% first tier
disp('Start constructing first tier...');
C_NTM_1 = 0.8;
[PPSC_x,PPSC_y] = ind2sub([img_height img_width],PPSC_index);
tier1 = delaunayTriangulation(PPSC_x,PPSC_y);
e = edges(tier1); 
n_e = size(e,1);
outliers = [];
d_s = -30:1:30;
d_v = -6:0.5:6;
d_k = -0.4:0.05:0.4;
n_s = size(d_s,2);
n_v = size(d_v,2);
n_k = size(d_k,2);
% if isempty(gcp('nocreate'))
%     parpool open local 8
% end
for i = 1:n_e
    P_ref_x = PPSC_x(e(i,1));
    P_ref_y = PPSC_y(e(i,1));
    P_end_x = PPSC_x(e(i,2));
    P_end_y = PPSC_y(e(i,2));
    y = squeeze(In(P_end_x,P_end_y,:).*conj(In(P_ref_x,P_ref_y,:)));
    norm_y = norm(y,2);
    
    S = 2*d_s.'*bn/(wave_length*center_slant_range);
    V = 2*d_v.'*b_t/(wave_length*1000);
    K = 2*d_k.'*b_k/(wave_length*1000);
    S = repmat(S,n_v,n_k);
    S = reshape(S.',baseline_number,n_k,n_s,n_v);
    S = permute(S,[3 4 2 1]);
    V = repmat(V,n_s,n_k);
    V = reshape(V.',baseline_number,n_k,n_v,n_s);
    V = permute(V,[4 3 2 1]);
    K = repmat(K,n_s,n_v);
    K = reshape(K.',baseline_number,n_v,n_k,n_s);
    K = permute(K,[4 2 3 1]);
    M = exp(-j*2*pi*(S+V+K));
    y = reshape(y,1,1,1,baseline_number);
    Gma = abs(sum(bsxfun(@times,M,y),4))/(sqrt(baseline_number)*norm_y);
    
%     for m = 1:size(d_s,2)
%         for n = 1:size(d_v,2)
%             for k = 1:size(d_k,2)
%                 a = exp(-j*2*pi*(2*bn*d_s(m)/(lamda*Rc)+2*b_t*d_v(n)/(lamda*1000)+2*b_k*d_k(k)/(lamda*1000)));
%                 Gma(m,n,k) = abs(a*y)/(sqrt(baseline_number)*norm_y);
%             end
%         end
%     end
    
    Gma_localmax = imregionalmax(Gma);
    Gma = Gma_localmax.*Gma;
    max_Gma = max(Gma(:));
    ind = find(Gma==max_Gma);
    if max_Gma < C_NTM_1
        outliers = [outliers,i];
    end
    if max_Gma >= C_NTM_1
        [e(i,3),e(i,4),e(i,5)] = ind2sub([n_s n_v n_k],ind(1));
        e(i,6) = max_Gma;
    end
end
% delete(gcp('nocreate'));
beep;
%--------------------------
fid = fopen('E:\data\TomoSAR\e.bin','w');
fwrite(fid,e,'double');
fclose(fid);
% 
% fid1 = fopen('E:\data\TomoSAR\e.bin');
% e = fread(fid1,[142 6],'*double');
% fclose(fid1);
%--------------------------
e(outliers,:)=[];
n_e = size(e,1);
e = findMaxSubnet(e);
%% M-estimator
disp('M-estimate for first tier...');
D=zeros(baseline_number,3);
D(:,1)=2*pi*2*bn/(wave_length*center_slant_range);
D(:,2)=2*pi*2*b_t/(wave_length*1000);
D(:,3)=2*pi*2*b_k/(wave_length*1000);
for i=1:n_e
    tic
    J=[d_s(e(i,3)),d_v(e(i,4)),d_k(e(i,5))];
    J=J.';
    delta_phi_1=D*J;
    P_ref_x = PPSC_x(e(i,1));
    P_ref_y = PPSC_y(e(i,1));
    P_end_x = PPSC_x(e(i,2));
    P_end_y = PPSC_y(e(i,2));
    phi_real = angle(squeeze(In(P_end_x,P_end_y,:).*conj(In(P_ref_x,P_ref_y,:))));
    delta_phi = delta_phi_1;
    for k=1:baseline_number
        wrapped = mod(delta_phi_1(k),2*pi);
        if abs(wrapped-phi_real(k))<pi
            delta_phi(k) = delta_phi(k)+phi_real(k)-wrapped;
        else if wrapped>phi_real(k)
                delta_phi(k) = delta_phi(k)+2*pi-(wrapped-phi_real(k));
            else
                delta_phi(k) = delta_phi(k)-(2*pi-(phi_real(k)-wrapped));
            end
        end
    end
    [e(i,3),e(i,4),e(i,5)] = M_estimate(D,delta_phi);
    toc
end
beep;
SPS_set = sort(unique([e(:,1);e(:,2)]));
SPS_index = PPSC_index(SPS_set);

%% network ajustment
disp('network ajustment for first tier...');
G = zeros(n_e,size(PPSC_index,1));
for i=1:n_e
    G(i,e(i,1))=-1;
    G(i,e(i,2))=1;
end
G = G(:,SPS_set);
G(:,1) = [];
I = eye(size(SPS_set,1)-1);
W = diag(e(:,6));
sigma=0:0.1:2;
E = zeros(1,size(sigma,2));
F = E;
SPS_1_list = [SPS_index,zeros(size(SPS_index,1),3)];
for i=3:5
    H = e(:,i);
    for k=1:size(sigma,2)
        X = (G.'*W*G+sigma(k)*I)\G.'*W*H;
        E(k) = norm(X);
        F(k) = norm(G*X-H);
    end
    kk = 1:1:size(sigma,2);
    figure;scatter3(F,E,kk);view(2);
    man_selected = input('Input the index of the corner');
    Sigma = sigma(man_selected);
    X = (G.'*W*G+Sigma*I)\G.'*W*H;
    SPS_1_list(2:size(SPS_1_list,1),i-1) = X;
end

%% second tier
disp('Start constucting second tier...');
window = 5;
C_NTM_2 = 0.75;
C_NTM_3 = 0.7;%0.6;
C_ratio = 0.8;%0.7;
[SPS_x,SPS_y] = ind2sub([img_height img_width],SPS_index);
dist = zeros(size(SPS_x,1),1);
SPSC_index = setdiff(SPSC_index,SPS_index);
SPS_2_index = [];
DPS_index = [];
count = 1;
for i=1:size(SPSC_index,1)
    [cur_x,cur_y] = ind2sub([img_height img_width],SPSC_index(i));
    for l=1:size(SPS_x,1)
        if(SPS_x(l)<=cur_x+window && SPS_x(l)>=cur_x-window && SPS_y(l)<=cur_y+window && SPS_y(l)>=cur_y-window)
            dist(l) = sqrt((cur_x-SPS_x(l))^2+(cur_y-SPS_y(l))^2);
        else
            dist(l) = Inf;
        end
    end
    if(min(dist)==Inf)
        
    else
        nearest = find(dist == min(dist));
        nearest_x = SPS_x(nearest);
        nearest_y = SPS_y(nearest);
        
        P_ref_x = nearest_x(1);
        P_ref_y = nearest_y(1);
        P_end_x = cur_x;
        P_end_y = cur_y;
        ref_index = sub2ind([img_height img_width],P_ref_x,P_ref_y);
        ref_params = SPS_1_list(SPS_1_list(:,1)==ref_index,2:4);
        y = squeeze(In(P_end_x,P_end_y,:).*conj(In(P_ref_x,P_ref_y,:)));
        norm_y = norm(y,2);
        S = 2*d_s.'*bn/(wave_length*center_slant_range);
        V = 2*d_v.'*b_t/(wave_length*1000);
        K = 2*d_k.'*b_k/(wave_length*1000);
        S = repmat(S,n_v,n_k);
        S = reshape(S.',baseline_number,n_k,n_s,n_v);
        S = permute(S,[3 4 2 1]);
        V = repmat(V,n_s,n_k);
        V = reshape(V.',baseline_number,n_k,n_v,n_s);
        V = permute(V,[4 3 2 1]);
        K = repmat(K,n_s,n_v);
        K = reshape(K.',baseline_number,n_v,n_k,n_s);
        K = permute(K,[4 2 3 1]);
        M = exp(-j*2*pi*(S+V+K));
        y = reshape(y,1,1,1,baseline_number);
        Gma = abs(sum(bsxfun(@times,M,y),4))/(sqrt(baseline_number)*norm_y);
%         for m = 1:size(d_s,2)
%             for n = 1:size(d_v,2)
%                 for k = 1:size(d_k,2)
%                     a = exp(j*2*pi*(2*bn*d_s(m)/(lamda*Rc)+2*b_t*d_v(n)/(lamda*1000)+2*b_k*d_k(k)/(lamda*1000)));
%                     Gma(m,n,k) = abs(conj(a)*y)/(sqrt(baseline_number)*norm(y,2));
%                 end
%             end
%         end
        Gma_localmax = imregionalmax(Gma);
        Gma = Gma_localmax.*Gma;
        max_Gma = max(Gma(:));
        ind = find(Gma==max_Gma);
        if max_Gma > C_NTM_2
            SPS_2_index = [SPS_2_index,i];
            [SPSC_index(i,2),SPSC_index(i,3),SPSC_index(i,4)] = ind2sub([n_s n_v n_k],ind(1));
            % M-estimate
            J=[d_s(SPSC_index(i,2)),d_v(SPSC_index(i,3)),d_k(SPSC_index(i,4))];
            J=J.';
            delta_phi_1=D*J;

            phi_real = angle(squeeze(In(P_end_x,P_end_y,:).*conj(In(P_ref_x,P_ref_y,:))));
            delta_phi = delta_phi_1;
            for o=1:baseline_number
                wrapped = mod(delta_phi_1(o),2*pi);
                if abs(wrapped-phi_real(o))<pi
                    delta_phi(o) = delta_phi(o)+phi_real(o)-wrapped;
                else if wrapped>phi_real(o)
                        delta_phi(o) = delta_phi(o)+2*pi-(wrapped-phi_real(o));
                    else
                        delta_phi(o) = delta_phi(o)-(2*pi-(phi_real(o)-wrapped));
                    end
                end
            end
            [SPSC_index(i,2),SPSC_index(i,3),SPSC_index(i,4)] = M_estimate(D,delta_phi);
            SPSC_index(i,2:4) = SPSC_index(i,2:4)+ref_params;
            %check for DPSs
            else if max_Gma > C_NTM_3
                    Gma(ind) = 0;
                    max_Gma2 = max(Gma(:));
                    if max_Gma2/max_Gma > C_ratio
                        ind2 = find(Gma==max_Gma2);
                        [SPSC_index(i,2),SPSC_index(i,3),SPSC_index(i,4)] = ind2sub([n_s n_v n_k],ind(1));
                        [SPSC_index(i,5),SPSC_index(i,6),SPSC_index(i,7)] = ind2sub([n_s n_v n_k],ind2(1));
                        SPSC_index(i,2:4) = [d_s(SPSC_index(i,2)),d_v(SPSC_index(i,3)),d_k(SPSC_index(i,4))];
                        SPSC_index(i,5:7) = [d_s(SPSC_index(i,5)),d_v(SPSC_index(i,6)),d_k(SPSC_index(i,7))];
                        if abs(SPSC_index(i,2)-SPSC_index(i,5))<50
                            DPS_index = [DPS_index,i];
                        end
                        SPSC_index(i,2:4) = SPSC_index(i,2:4)+ref_params;
                        SPSC_index(i,5:7) = SPSC_index(i,5:7)+ref_params;
                    end
                end
        end
    end
    disp(count);
    count = count+1;
end
beep;
if size(SPSC_index,2)==4
    SPSC_index(1,5:7) = [0 0 0];
end
SPS_2_list = SPSC_index(SPS_2_index,1:4);
DPS_list = SPSC_index(DPS_index,1:7);
%--------------------------
% fid = fopen('E:\data\TomoSAR\SPS_2_list_07_08.bin','w');
% fwrite(fid,SPS_2_list,'double');
% fclose(fid);
% 
% fid1 = fopen('E:\data\TomoSAR\SPS_2_list_07_08.bin');
% SPS_2_list = fread(fid1,[182 4],'*double');
% fclose(fid1);
%--------------------------
%--------------------------
% fid = fopen('E:\data\TomoSAR\DPS_list_07_08.bin','w');
% fwrite(fid,DPS_list,'double');
% fclose(fid);
% 
% fid1 = fopen('E:\data\TomoSAR\DPS_list_07_08.bin');
% DPS_list = fread(fid1,[39 7],'*double');
% fclose(fid1);
%--------------------------
%% Display
if size(SPS_1_list,1)==0
    x1 =[]; y1 = [];vertical_h1 = [];v1 = [];k1 = [];
else
    [x1,y1] = ind2sub([img_height img_width],SPS_1_list(:,1));
    vertical_h1 = SPS_1_list(:,2)*sin(incidence_angle);
    x1 = x1*azimuth_spacing;
    y1 = (y1+SPS_1_list(:,2)*sin(incidence_angle)*cos(incidence_angle)/range_spacing)*range_spacing/sin(incidence_angle);
    v1 = SPS_1_list(:,3);
    k1 = SPS_1_list(:,4);
end
if size(SPS_2_list,1)==0
    x2 = [];y2 = [];vertical_h2 = [];v2 = [];k2 = [];
else
    [x2,y2] = ind2sub([img_height img_width],SPS_2_list(:,1));
    vertical_h2 = SPS_2_list(:,2)*sin(incidence_angle);
    x2 = x2*azimuth_spacing;
    y2 = (y2+SPS_2_list(:,2)*sin(incidence_angle)*cos(incidence_angle)/range_spacing)*range_spacing/sin(incidence_angle);
    v2 = SPS_2_list(:,3);
    k2 = SPS_2_list(:,4);
end
if size(DPS_index,1) == 0
    x3 = [];y3 = [];vertical_h3 = [];v3 = [];k3 = [];
    x4 = [];y4 = [];vertical_h4 = [];v4 = [];k4 = [];
else
    [x3,y3] = ind2sub([img_height img_width],DPS_list(:,1));
    x3 = x3*azimuth_spacing;
    x4 = x3;
    vertical_h3 = DPS_list(:,2)*sin(incidence_angle);
    vertical_h4 = DPS_list(:,5)*sin(incidence_angle);
    y4 = (y3+DPS_list(:,5)*sin(incidence_angle)*cos(incidence_angle)/range_spacing)*range_spacing/sin(incidence_angle);
    y3 = (y3+DPS_list(:,2)*sin(incidence_angle)*cos(incidence_angle)/range_spacing)*range_spacing/sin(incidence_angle);
    v3 = DPS_list(:,3);k3 = DPS_list(:,6);
    v4 = DPS_list(:,4);k4 = DPS_list(:,7);
end
%height
figure;
scatter3(x1,y1,vertical_h1,10,vertical_h1,'+');
hold on
scatter3(x2,y2,vertical_h2,10,vertical_h2,'d');
hold on
scatter3(x3,y3,vertical_h3,5,vertical_h3,'s');
hold on
scatter3(x4,y4,vertical_h4,10,vertical_h4,'x');
%deformation
figure;
scatter3(x1,y1,vertical_h1,10,v1,'+');
hold on
scatter3(x2,y2,vertical_h2,10,v2,'d');
hold on
scatter3(x3,y3,vertical_h3,5,v3,'s');
hold on
scatter3(x4,y4,vertical_h4,10,v4,'x');
%temperature
figure;
scatter3(x1,y1,vertical_h1,10,k1,'+');
hold on
scatter3(x2,y2,vertical_h2,10,k2,'d');
hold on
scatter3(x3,y3,vertical_h3,5,k3,'s');
hold on
scatter3(x4,y4,vertical_h4,10,k4,'x');
