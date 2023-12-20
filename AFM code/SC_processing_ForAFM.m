[FileName,PathName] = uigetfile('*.tsv'); 
Table_raw = readtable([PathName,FileName],"FileType","text",'PreserveVariableNames',true); %variables text

[FileName,PathName] = uigetfile('*.xlsx'); 
Table_raw = readtable([PathName,FileName]);

%% 筛选raw列表中要用到的列,存入新表Table_selected,重新赋变量名称 filtering column from table_raw and store it into Table_selected, rename the varioable names of column.
Table_selected = Table_raw(:,[1;8;16;23;29;30;31;38;44;45;46;53;59;60;61;68;74;75;76;83;89;90;91;98;104;105;106]);
Table_selected.Properties.VariableNames = {'Filename','YM','Height',...
    'IM_1hz','G_prime_1hz','G_doubleprime_1hz','LT_1hz',...
    'IM_3hz','G_prime_3hz','G_doubleprime_3hz','LT_3hz',...
    'IM_10hz','G_prime_10hz','G_doubleprime_10hz','LT_10hz',...
    'IM_30hz','G_prime_30hz','G_doubleprime_30hz','LT_30hz',...
    'IM_100hz','G_prime_100hz','G_doubleprime_100hz','LT_100hz',...
    'IM_200hz','G_prime_200hz','G_doubleprime_200hz','LT_200hz'};
%% 决定测量每个细胞时的曲线数量，并加标签（需确认每个细胞的坐标改变,每个细胞重复测量的坐标不变）decide for how many curves are measured for each single cells, and put an index on it(make sure the coordinate is changed)
X_Position = Table_raw.('X Position');
Y_Position = Table_raw.('Y Position');
Cell_position_x = X_Position(1);
Cell_position_y = Y_Position(1);
sz = size(Table_selected,1);
Index_m = zeros(sz,2);
j=1;
k=1;
for i=1:sz
    Cell_position_xx= X_Position(i);
    Cell_position_yy= Y_Position(i);
    if Cell_position_xx == Cell_position_x && Cell_position_yy == Cell_position_y
        Index_m(i,1) = k;
        Index_m(i,2) = j;
        j=j+1;
    else
        Cell_position_x = X_Position(i);
        Cell_position_y = Y_Position(i);
        k=k+1;
        Index_m(i,1) = k;
        Index_m(i,2) = 1;
        j=2;
    end
end
Table_selected.Index_percell = Index_m(:,1);
Table_selected.Index_percurve = Index_m(:,2);
Table_selected = movevars(Table_selected, 'Index_percurve', 'Before', 'Filename');
Table_selected = movevars(Table_selected, 'Index_percell', 'Before', 'Index_percurve');
clearvars  -except PathName FileName Table_raw Table_selected Matrix_percell_Results Table_Final_Results

%% 把table转换为矩阵，去除第一列文件名，存为Matrix_cellmech     table to matrix, delete the first name row, save it into Matrix_cellmech
% Matrix_cellmech = table2array(Cell_mech(:,2:27));

%% Filtering
sz = size(Table_selected,1); %总数据sz
a = Table_selected.Index_percell(sz); %每五条曲线为一个细胞，取a个细胞
Matrix_percell_Results = zeros(a,17); %创建写入结果的矩阵，第一列为编号，第二列为每个细胞的重复测量数，第三列为每个细胞的第一条曲线编号，第四列开始自定（此处为中位数）
Matrix_percell_Results(:,1) = [1:a]; %给每个细胞编号 put an index for each cell
I = find(Table_selected.Index_percurve == 1); 
Matrix_percell_Results(:,3) = I;
for i = 1:(a-1)
     I(i) = I(i+1) - I(i);
end
I(a) = sz - I(a) + 1;
Matrix_percell_Results(:,2) = I;  %结果矩阵第二列输入每个细胞的重复测量数 the 2nd column is the number of curves for repeating measure on one cell
Matrix_percell_Results(:,3) = Matrix_percell_Results(:,3) + I; %the 3rd column is the index of first curve on cell
clearvars  -except PathName FileName Table_raw Table_selected Matrix_percell_Results Table_Final_Results sz a
%% YM的中位数  median of YM
j=1;
k=1;
Matrix_percell_YM = zeros(Matrix_percell_Results(1,2),1);
for i=1:sz
    if i ~= Matrix_percell_Results(j,3) && i ~= sz
        Matrix_percell_YM(k,1) = Table_selected.YM(i);  %此处循环让Matrix_percell_YM暂时存储每个细胞的对应数据（此处为YM）
        k=k+1; 
    elseif i == sz
        Matrix_percell_YM(k,1) = Table_selected.YM(i);
        Matrix_percell_Results(j,4) = median(Matrix_percell_YM);
    else
        Matrix_percell_Results(j,4) = median(Matrix_percell_YM); %矩阵Matrix_percell第三列为每个细胞的YM中位数
        j=j+1;
        Matrix_percell_YM = zeros(Matrix_percell_Results(j,2),1);
        k=2;
        Matrix_percell_YM(1,1) = Table_selected.YM(i);
    end 
end
clearvars -except FileName PathName Table_raw Table_selected Matrix_percell_Results Table_Final_Results sz a
%% 挑选G'和G''的曲线

for i=1:a
    Table_percell=Table_selected(Table_selected.Index_percell==i,:);
    [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
    Matrix_G_complex = [Table_percell.G_prime_1hz(Median_idx) Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx) Table_percell.G_prime_200hz(Median_idx);...
                        Table_percell.G_doubleprime_1hz(Median_idx) Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx) Table_percell.G_doubleprime_200hz(Median_idx)];
    TF = issorted(Matrix_G_complex,2); %检验是否为升序 if ascending 
    count = 1;
    while TF == false
       Table_percell(Median_idx,:) = [];
       [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
       Matrix_G_complex = [Table_percell.G_prime_1hz(Median_idx) Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx) Table_percell.G_prime_200hz(Median_idx);...
                          Table_percell.G_doubleprime_1hz(Median_idx) Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx) Table_percell.G_doubleprime_200hz(Median_idx)];
       TF = issorted(Matrix_G_complex,2);
       count = count + 1;
    end                    
       if count <= Matrix_percell_Results(i,2) 
            Matrix_percell_Results(i,5) = Dif_Median;
            Matrix_percell_Results(i,[6:11]) = Matrix_G_complex(1,[1:6]); 
            Matrix_percell_Results(i,12) = Table_selected.Index_percurve(Median_idx);
            Matrix_percell_Results(i,[13:18]) = Matrix_G_complex(2,[1:6]);
       else   %如果都不符合升序的话，排除1hz的数据 if not ascending, delete 1hz
               Table_percell=Table_selected(Table_selected.Index_percell==i,:);
               [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
               Matrix_G_complex = [Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx) Table_percell.G_prime_200hz(Median_idx);...
                                   Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx) Table_percell.G_doubleprime_200hz(Median_idx)];
               TF = issorted(Matrix_G_complex,2);
               count = 1;
               while TF == false
                   Table_percell(Median_idx,:) = [];
                   [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
                   Matrix_G_complex = [Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx) Table_percell.G_prime_200hz(Median_idx);...
                                       Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx) Table_percell.G_doubleprime_200hz(Median_idx)];
                   TF = issorted(Matrix_G_complex,2);
                   count = count + 1;
               end
               if count <= Matrix_percell_Results(i,2)
                   Matrix_percell_Results(i,5) = Dif_Median;
                   Matrix_percell_Results(i,6) = 0;
                   Matrix_percell_Results(i,[7:11]) = Matrix_G_complex(1,[1:5]); 
                   Matrix_percell_Results(i,12) = Table_selected.Index_percurve(Median_idx);
                   Matrix_percell_Results(i,13) = 0;
                   Matrix_percell_Results(i,[14:18]) = Matrix_G_complex(2,[1:5]);
               else   %如果排除1hz的数据都不符合升序的话，排除200hz if not ascending without 1hz, delete 200hz
                       Table_percell=Table_selected(Table_selected.Index_percell==i,:);
                       [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
                       Matrix_G_complex = [Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx);...
                                           Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx)];
                       TF = issorted(Matrix_G_complex,2);
                       count = 1;
                   while TF == false
                       Table_percell(Median_idx,:) = [];
                       [Dif_Median, Median_idx] = min(abs(Table_percell.YM - Matrix_percell_Results(i,4)));
                       Matrix_G_complex = [Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx);...
                                           Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx)];
                       TF = issorted(Matrix_G_complex,2);
                       count = count + 1;
                   end
                   if count <= Matrix_percell_Results(i,2)
                       Matrix_percell_Results(i,5) = Dif_Median;
                       Matrix_percell_Results(i,6) = 0;
                       Matrix_percell_Results(i,[7:10]) = Matrix_G_complex(1,[1:4]); 
                       Matrix_percell_Results(i,11) = 0;
                       Matrix_percell_Results(i,12) = Table_selected.Index_percurve(Median_idx);
                       Matrix_percell_Results(i,13) = 0;
                       Matrix_percell_Results(i,[14:17]) = Matrix_G_complex(2,[1:4]);
                       Matrix_percell_Results(i,18) = 0;
               else   %如果排除1hz的数据都不符合升序的话，剔除数据 if still not working, delete the whole data
                   end
               end
       end
end
clearvars -except FileName PathName Table_raw Table_selected Matrix_percell_Results Table_Final_Results
Table_Final_Results = array2table(Matrix_percell_Results);
Table_Final_Results.Properties.VariableNames = {'Cell_Index','Curves_percell','Next_Cell_Index','YM','Dif_YM'...
    'G_prime_1hz','G_prime_3hz','G_prime_10hz','G_prime_30hz','G_prime_100hz','G_prime_200hz','Median_Index',...
    'G_doubleprime_1hz','G_doubleprime_3hz','G_doubleprime_10hz','G_doubleprime_30hz','G_doubleprime_100hz','G_doubleprime_200hz'};

%% 临时筛选 manual filtering 
    Table_percell=Table_selected(Table_selected.Index_percell==43,:); %改index数字 change index number manually
    Matrix_G_complex_total=[];
for Median_idx=1:size(Table_percell,1)
    Matrix_G_complex = [Table_percell.G_prime_1hz(Median_idx) Table_percell.G_prime_3hz(Median_idx) Table_percell.G_prime_10hz(Median_idx) Table_percell.G_prime_30hz(Median_idx) Table_percell.G_prime_100hz(Median_idx) Table_percell.G_prime_200hz(Median_idx);...
                        Table_percell.G_doubleprime_1hz(Median_idx) Table_percell.G_doubleprime_3hz(Median_idx) Table_percell.G_doubleprime_10hz(Median_idx) Table_percell.G_doubleprime_30hz(Median_idx) Table_percell.G_doubleprime_100hz(Median_idx) Table_percell.G_doubleprime_200hz(Median_idx)];
    TF = issorted(Matrix_G_complex,2);
    Matrix_G_complex_total = [Matrix_G_complex_total;Matrix_G_complex];
end
%% FE model fitting
data = table2struct(Table_raw);

% for k=1:length(data)
%     data2(k).f = [1 3 10 30 100 200];
%     data2(k).G1 = [data(k).G_prime_1hz data(k).G_prime_3hz data(k).G_prime_10hz data(k).G_prime_30hz data(k).G_prime_100hz data(k).G_prime_200hz];
%     data2(k).G2 = [data(k).G_doubleprime_1hz data(k).G_doubleprime_3hz data(k).G_doubleprime_10hz data(k).G_doubleprime_30hz data(k).G_doubleprime_100hz data(k).G_doubleprime_200hz];
% end


%beads data from Cath G'=G__   G''=G___
% for k=1:length(data)
%     data2(k).f = [1 3 10 30 100 200];
%     data2(k).G1 = [data(k).G__1Hz data(k).G__3Hz data(k).G__10Hz data(k).G__30Hz data(k).G__100Hz data(k).G__200Hz];
%     data2(k).G2 = [data(k).G___1Hz data(k).G___3Hz data(k).G___10Hz data(k).G___30Hz data(k).G___100Hz data(k).G___200Hz];
% end
for k=1:length(data)
    data2(k).f = [1 3 10 30 100];
    data2(k).G1 = [data(k).G__1Hz data(k).G__3Hz data(k).G__10Hz data(k).G__30Hz data(k).G__100Hz];
    data2(k).G2 = [data(k).G___1Hz data(k).G___3Hz data(k).G___10Hz data(k).G___30Hz data(k).G___100Hz];
end


% for k=1:length(data)
%     data2(k).f = [1 3 10 30 100];
%     data2(k).G1 = [data(k).G_prime_1hz data(k).G_prime_3hz data(k).G_prime_10hz data(k).G_prime_30hz data(k).G_prime_100hz];
%     data2(k).G2 = [data(k).G_doubleprime_1hz data(k).G_doubleprime_3hz data(k).G_doubleprime_10hz data(k).G_doubleprime_30hz data(k).G_doubleprime_100hz];
% end

MasterList = data2;

for k= 1:length (MasterList)
    [mu, alpha]= springpotfit(MasterList(k).f,MasterList(k).G1,MasterList(k).G2);
    MasterList(k).mu = mu;
    MasterList(k).alpha = alpha;
end


% load MasterList %maybe add your name
%select frequency for G1 and G2 selcetion only
Table_selected = Table_Final_Results_10A; 
Table_selected = Table_Final_Results_436;

MasterList_temp = table2struct(Table_selected);

for i = 1:length(MasterList_temp)
MasterList(i).mu = MasterList_temp(i).mu_1to100;
MasterList(i).alpha = MasterList_temp(i).alpha_1to100;
end

t= 0.1; %select a frequency for calculation of G*  (and G' and G'') from mu and alpha...dont know why it is called t here...f was taken
% Gstarc = Magnitude or absolute of G* calculated from mu and alpha 
clear i

for k = 1:length(MasterList)
                    MasterList(k).Gstarc=abs(MasterList(k).mu^(1-MasterList(k).alpha)*(i*2*pi*t)^MasterList(k).alpha); % This is your Stiffness or Resistance based on mu and alpha from the whole frequency range
                    MasterList(k).phi = MasterList(k).alpha * pi/2; %Fluidity
                    MasterList(k).G1c = real(MasterList(k).mu^(1-MasterList(k).alpha)*(i*2*pi*t)^MasterList(k).alpha); 
                    MasterList(k).G2c = imag(MasterList(k).mu^(1-MasterList(k).alpha)*(i*2*pi*t)^MasterList(k).alpha);  
end

MasterList_table = struct2table(MasterList);
MasterList_table = MasterList_table (:,[3:6]);
Table_selected =cat(2,Table_selected, MasterList_table);

Table_Final_Results_10A = Table_selected;
save(['D:\Documents\AFM data\SC-AFM-STAINED\Results','/Table_Final_Results_10A_30hzGstar.mat'],'Table_Final_Results_10A')
Table_Final_Results_436 = Table_selected;
save(['D:\Documents\AFM data\SC-AFM-STAINED\Results','/Table_Final_Results_436_30hzGstar.mat'],'Table_Final_Results_436')
clearvars -except Table_Final_Results_10A Table_Final_Results_436

clearvars -except FileName PathName Table_raw Table_selected Matrix_percell_Results Table_Final_Results
%% 导入AR,Area add AR.Area
[FileName_Shape,PathName_Shape] = uigetfile('*.csv'); 
Table_Shape_raw = readtable([PathName_Shape,FileName_Shape],"FileType","text",'PreserveVariableNames',true,'ReadVariableNames',true); %变量为文本
Table_Shape_raw = Table_Shape_raw(:,{'Area','AR'});
Table_Final_Results = [Table_Final_Results Table_Shape_raw];
clearvars -except FileName PathName Table_Final_Results


%% 保存结果10A save results

Table_FR_Selected_10A = Table_Final_Results(:,{'Cell_Index','Area','AR','YM','mu','alpha'});
save([PathName,'/Table_Final_Results_10A.mat'],'Table_Final_Results')
save([PathName,'/Table_FR_Selected_10A.mat'],'Table_FR_Selected_10A')
save(['D:\Documents\Viscoelasticity Measurment\SC-AFM-STAINED\Results','/Table_FR_Selected_10A_fixed.mat'],'Table_Final_Results_10A')

%% 保存结果436 save results

Table_FR_Selected_436 = Table_Final_Results(:,{'Cell_Index','Area','AR','YM','mu','alpha'});
save([PathName,'/Table_Final_Results_436.mat'],'Table_Final_Results')
save([PathName,'/Table_FR_Selected_436.mat'],'Table_FR_Selected_436')
save(['D:\Documents\Viscoelasticity Measurment\SC-AFM-STAINED\Results','/Table_FR_Selected_436_fixed.mat'],'Table_Final_Results_436')

%% plot，distribution
distributionFitter
x=Table_FR_Selected_10A.AR;
y=Table_FR_Selected_10A.mu;
y1=Table_FR_Selected_10A.alpha;
x1=Table_FR_Selected_436.AR;
y2=Table_FR_Selected_436.mu;
y3=Table_FR_Selected_436.alpha;
labs=num2str(data_to_delete.Var1);

figure
plot(x,y, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
hold on
plot(x1,y1, 'b+', 'MarkerSize', 15, 'LineWidth', 2);

figure
plot(x,y1, 'ro', 'MarkerSize', 15, 'LineWidth', 2);
hold on
plot(x1,y3, 'bo', 'MarkerSize', 15, 'LineWidth', 2);
%% 平均值 average
mean_mu_10A_lowAR = mean(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR < 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));
median_mu_10A_lowAR = median(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR < 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));
std_mu_10A_lowAR = std(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR < 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));

mean_mu_10A_highAR = mean(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR > 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));
median_mu_10A_highAR = median(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR > 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));
std_mu_10A_highAR = std(Table_Final_Results_10A.mu_1to100(Table_Final_Results_10A.AR > 1.2 & Table_Final_Results_10A.mu_1to100 ~= 0 ,:));