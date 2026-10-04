clear, close all
MX_Switch_SG = 'MX';
tmotload = '4';
isotope = '162';
n_count = 0;
MX_Switch_SG_cell = {'MX';'Switch';'SG'};
n_MX_Switch_SG_cell = length(MX_Switch_SG_cell);
for k_MX_Switch_SG_cell = 1:n_MX_Switch_SG_cell
    MX_Switch_SG = MX_Switch_SG_cell{k_MX_Switch_SG_cell};
    for tmotload_num = [4 8 12 16 20 24 30]%['4';'8';'12';'16';'20';'24';'30']
        tmotload = num2str(tmotload_num);
        for isotope_num = [162 164]
            isotope = num2str(isotope_num);
            n_count = n_count +1;
path1 = ['\\10.16.18.204\实验室内部文件\Results\MOT_single_vs_mixture\20260726处理0712\Temperature\',MX_Switch_SG,'\',isotope,'\tmotload',tmotload,' ',MX_Switch_SG,isotope];
load([path1,'\Mat.mat'],'sigmax','sigmay','toftime','sigmax_std','sigmay_std')
toftime_ms = toftime*1e3;
weight_simgax = 1./sigmax_std;
weight_simgay = 1./sigmay_std;

k_B = 1.380649e-23;
amu = 1.66053906660e-27;
m162 = 161.926795;
m164 = 163.929171;
% --------------------x-------------------------------------------------------------------------------------------------------------------------
[fitresultx, gofx] = cloudTfit(toftime, sigmax, weight_simgax,eval(['m',isotope]),'x',isotope,tmotload,MX_Switch_SG);
figure(888)
errorbar(toftime_ms,sigmax*1e6,sigmax_std*1e6,sigmax_std*1e6,...
    'LineStyle','none','CapSize',0,'Marker','o',...
    'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',3);
hold on
plot(toftime_ms,fitresultx(toftime)*1e6,'LineWidth',1,'Color','r')
title([MX_Switch_SG,', ',isotope,', load:',tmotload,'s'])
xlabel('Time of flight (ms)')
ylabel('\sigma_x (um)')
LGD = legend('Raw data','fitting');
LGD.Box = "off";
LGD.Position = [0.6 0.2 0.3 0.15];
% ---- 设置字体为 Times New Roman ----
ax = gca;                           
ax.FontName = 'Times New Roman';    % 刻度标签字体
ax.FontSize = 7;                    % ★ 刻度标签字体大小为 7
ax.XLabel.FontName = 'Times New Roman';  
ax.XLabel.FontSize = 8;             % ★ X轴标题字体大小为 8
ax.YLabel.FontName = 'Times New Roman';  
ax.YLabel.FontSize = 8;             % ★ Y轴标题字体大小为 8
% ---- 设置图形尺寸 ----
fig = gcf;
fig.Units = 'centimeters';
fig.Position(3:4) = [8.5, 4.3];
% ---- 导出为矢量 PDF ----
exportgraphics(fig, [path1,'\TOF_fitx(new).pdf'], 'ContentType', 'vector');
saveas(fig, [path1,'\TOF_fitx(new).svg']);
exportgraphics(fig, [path1,'\TOF_fitx(new).png'], 'Resolution', 600);
% --------------------y-------------------------------------------------------------------------------------------------------------------------
[fitresulty, gofy] = cloudTfit(toftime, sigmay, weight_simgay,eval(['m',isotope]),'y',isotope,tmotload,MX_Switch_SG);
figure(999)
errorbar(toftime_ms,sigmay*1e6,sigmay_std*1e6,sigmay_std*1e6,...
    'LineStyle','none','CapSize',0,'Marker','o',...
    'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',3);
hold on
plot(toftime_ms,fitresulty(toftime)*1e6,'LineWidth',1,'Color','b')
title([MX_Switch_SG,', ',isotope,', load:',tmotload,'s'])
xlabel('Time of flight (ms)')
ylabel('\sigma_y (um)')
LGD = legend('Raw data','fitting');
LGD.Box = "off";
LGD.Position = [0.6 0.2 0.3 0.15];
% ---- 设置字体为 Times New Roman ----
ax = gca;                           
ax.FontName = 'Times New Roman';    % 刻度标签字体
ax.FontSize = 7;                    % ★ 刻度标签字体大小为 7
ax.XLabel.FontName = 'Times New Roman';  
ax.XLabel.FontSize = 8;             % ★ X轴标题字体大小为 8
ax.YLabel.FontName = 'Times New Roman';  
ax.YLabel.FontSize = 8;             % ★ Y轴标题字体大小为 8
% ---- 设置图形尺寸 ----
fig = gcf;
fig.Units = 'centimeters';
fig.Position(3:4) = [8.5, 4.3];
% ---- 导出为矢量 PDF ----
exportgraphics(fig, [path1,'\TOF_fity(new).pdf'], 'ContentType', 'vector');
saveas(fig, [path1,'\TOF_fity(new).svg']);
exportgraphics(fig, [path1,'\TOF_fity(new).png'], 'Resolution', 600);
%-------------------后处理-----------------------
cix = confint(fitresultx, 0.95);
ciy = confint(fitresulty, 0.95);
Tx = fitresultx.T;
Ty = fitresulty.T;
dfx = gofx.dfe;
t_val = tinv(0.975, dfx); % 对于95% CI，两尾
SE_Tx = (cix(2,1)-cix(1,1))/(2*t_val);
SE_Ty = (ciy(2,1)-ciy(1,1))/(2*t_val);

sigmax0 = fitresultx.sigmax0;
sigmay0 = fitresulty.sigmax0;
SE_sigmax0 = (cix(2,2) - cix(1,2)) / (2 * t_val);
SE_sigmay0 = (ciy(2,2) - ciy(1,2)) / (2 * t_val);

% ========== 写入 Tfit.txt ==========
fid = fopen([path1,'\Tfit(new).txt'], 'w');
fprintf(fid, '========== 温度测量结果 ==========\n');
fprintf(fid, 'Tx = %.2f ± %.2f µK \n', Tx, SE_Tx);
fprintf(fid, 'Ty = %.2f ± %.2f µK \n', Ty, SE_Ty);
fprintf(fid, '\n');
fprintf(fid, '========== 初始尺寸 (sigma0) ==========\n');
fprintf(fid, 'sigmax0 = %.3f ± %.3f µm \n', sigmax0*1e6, SE_sigmax0*1e6);
fprintf(fid, 'sigmay0 = %.3f ± %.3f µm \n', sigmay0*1e6, SE_sigmay0*1e6);
fprintf(fid, '=======================================\n');
fclose(fid);
% ========== save mat ==========
save([path1,'\Mat_new.mat'],'fitresultx','fitresulty','gofx','gofy','Tx','Ty','sigmax0','sigmax0',...
    'SE_Tx','SE_Ty','SE_sigmax0','SE_sigmay0',...
    'toftime','sigmax','sigmay','sigmax_std','sigmay_std',...
    'm162','m164')
close all
        end
    end
end

% ---------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------------------

function [fitresult, gof] = cloudTfit(toftime, sigmax, weight_simgax,mass,xory,isotope,tmotload,MX_Switch_SG)

[xData, yData, weights] = prepareCurveData( toftime, sigmax, weight_simgax );

% set fittype & options
fittype_text = ['sqrt(sigmax0^2+1.380649e-2/(',num2str(mass,'%.6f'),'*1.66053906660)*T*t.^2)'];
ft = fittype( fittype_text, 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
if xory == 'x'
    opts.StartPoint = [7 0.032];
else
    opts.StartPoint = [9 0.015];
end
opts.Weights = weights;

% fitting
[fitresult, gof] = fit( xData, yData, ft, opts );
disp(['fitting finished: ',MX_Switch_SG,isotope,' ',tmotload,' sigma',xory])
disp(['sigma0',xory,': ',num2str(fitresult.sigmax0),' m'])
disp(['T: ',num2str(fitresult.T),' uK'])
disp(['R2: ',num2str(gof.rsquare)])
disp('--------------------')
end


