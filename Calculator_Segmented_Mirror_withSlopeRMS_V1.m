clear; clc; close all;

%% 1. 数据导入和预处理
data = readmatrix('OS_LOS0_Local.txt');

nodeID = data(:,1);
X = data(:,2);
Y = data(:,3);
Z = data(:,4);
UX = data(:,5);
UY = data(:,6);
UZ = data(:,7);

% 变形后节点
Xd = X + UX;
Yd = Y + UY;
Zd = Z + UZ;

orig_points = [X, Y, Z];
deformed_points = [Xd, Yd, Zd];

%% 2. 刚体运动最小二乘移除（Kabsch）
center_orig = mean(orig_points,1);
center_def = mean(deformed_points,1);
P1 = orig_points - center_orig;
P2 = deformed_points - center_def;

H = P1' * P2;
[U, ~, V] = svd(H);
R = V * U';
if det(R) < 0
    V(:,3) = -V(:,3);
    R = V * U';
end
T = center_def' - R * center_orig';
P_corr = (R * orig_points' + T)';    % 矫正后的刚体位置
W = deformed_points - P_corr;        % 残差（面型误差） Nx3

% 明确 Wz 为法向位移（第三列），以 m 为单位
% 按你的脚本原样：使用 UZ 作为 Wz
%Wz = W(:,3);
Wz = UZ;

[pitch,roll,yaw]=dcm2angle(R,"XYZ");

PV = (max(Wz) - min(Wz)) * 1e6;
RMS = rms(Wz) * 1e6;

RigidBody_Motions=[T(1)*1e3;T(2)*1e3;T(3)*1e3;pitch*1e6;roll*1e6;yaw*1e6];

RigidBody_Iterm={'Tx (um)';'Ty (um)';'Tz (um)';'Rx (urad)';'Ry (urad)';'Rz (urad)'};
RigidBody_Matrix=table(RigidBody_Iterm,RigidBody_Motions);
RigidBody_Matrix.RigidBody_Iterm=char(RigidBody_Matrix.RigidBody_Iterm);
fprintf('\n===光学表面刚体运动===\n\n');
disp(RigidBody_Matrix)

% fprintf('\n===刚体运动===\n');
% fprintf('Tx (um):%.6\n', T(1)*1e3);
% fprintf('Ty (um):%.6\n', T(2)*1e3);
% fprintf('Tz (um):%.6\n', T(3)*1e3);
% fprintf('Rx (urad):%.3\n', pitch*1e6);
% fprintf('Ry (urad):%.3\n', roll*1e6);
% fprintf('Rz (urad):%.3\n', yaw*1e6);
%fprintf('刚体运动移除后：PV=%.2f nm, RMS=%.2f nm\n', PV, RMS);

%% 3. Zernike多项式拟合
zernike_order = 10; % 可调整
rho = sqrt(X.^2 + Y.^2) / max(sqrt(X.^2 + Y.^2)); % 归一到单位圆
theta = atan2(Y, X);
Zmat = generateZernikeBasis(rho, theta, zernike_order); % N x M

[Zmat_ortho, ~] = qr(Zmat, 0); % 列正交基

zernike_coef = Zmat_ortho \ Wz;
zernike_fit = Zmat_ortho * zernike_coef;        % m (单位与 Wz 一致)
zernike_resid = Wz - zernike_fit;               % m

%fprintf('Zernike 拟合完成：系数个数 = %d\n', length(zernike_coef));
%fprintf('Zernike RMS残差（含所有模态）：%.2f nm\n', rms(zernike_resid) * 1e6);

% ---------- 去除 PTT 与 PTT+Power（Zernike） ----------
modes = [];
for n = 0:zernike_order
    for m = -n:2:n
        modes = [modes; n, m];
    end
end
ptt_idx = 1:min(3, size(Zmat_ortho,2)); % piston, tip, tilt 列索引（生成顺序一致）
defocus_idx = find(modes(:,1)==2 & modes(:,2)==0);
do_defocus_zern = ~isempty(defocus_idx);

% PTT term（来自 fit）
PTT_term_zern = Zmat_ortho(:, ptt_idx) * zernike_coef(ptt_idx); % 与 Wz 单位一致
resid_ptt_removed_zern = Wz - PTT_term_zern;                     % 与 Wz 单位一致

RMS_Zern_PTT_Removed = rms(resid_ptt_removed_zern);
PV_Zern_PTT_Removed = max(resid_ptt_removed_zern) - min(resid_ptt_removed_zern);
%fprintf('（Zernike）去除 PTT 后 RMS=%.2f nm, PV=%.2f nm\n', RMS_Zern_PTT_Removed*1e6, PV_Zern_PTT_Removed*1e6);

if do_defocus_zern
    idx_ptt_power = unique([ptt_idx(:); defocus_idx(:)]);
    PTTPower_term_zern = Zmat_ortho(:, idx_ptt_power) * zernike_coef(idx_ptt_power); % same unit
    resid_pttpower_removed_zern = Wz - PTTPower_term_zern;                            % same unit
    RMS_Zern_PTTPower_Removed = rms(resid_pttpower_removed_zern);
    PV_Zern_PTTPower_Removed = max(resid_pttpower_removed_zern) - min(resid_pttpower_removed_zern);
    %fprintf('（Zernike）去除 PTT+Power 后 RMS=%.2f nm, PV=%.2f nm\n', RMS_Zern_PTTPower_Removed*1e6, PV_Zern_PTTPower_Removed*1e6);
else
    PTTPower_term_zern = zeros(size(Wz));
    resid_pttpower_removed_zern = Wz;
end

%% 4. XY多项式拟合
xy_order = 10;
[XYmat, xy_terms] = generateXYPolynomialBasis(X, Y, xy_order); % N x M, xy_terms: [i j] per column

[XYmat_ortho, ~] = qr(XYmat, 0);

xy_coef = XYmat_ortho \ Wz;
xy_fit = XYmat_ortho * xy_coef;   % same unit as Wz
xy_resid = Wz - xy_fit;           % same unit as Wz

%fprintf('XY 多项式 拟合完成：系数个数 = %d\n', length(xy_coef));
%fprintf('XY 多项式 RMS残差（含所有项）：%.2f nm\n', rms(xy_resid) * 1e6);

% ---------- 去除 PTT 与 PTT+Power（XY） ----------
ptt_idx_xy = 1:min(3, size(XYmat_ortho,2));
PTT_term_xy = XYmat_ortho(:, ptt_idx_xy) * xy_coef(ptt_idx_xy); % same unit
resid_ptt_removed_xy = Wz - PTT_term_xy;                        % same unit

RMS_XY_PTT_Removed = rms(resid_ptt_removed_xy);
PV_XY_PTT_Removed = max(resid_ptt_removed_xy) - min(resid_ptt_removed_xy);
%fprintf('（XY）去除 PTT 后 RMS=%.2f nm, PV=%.2f nm\n', RMS_XY_PTT_Removed*1e6, PV_XY_PTT_Removed*1e6);

% 构造 power 向量并与 PTT 正交化
power_vec = (X.^2 + Y.^2);
p = power_vec(:);
proj_on_ptt = XYmat_ortho(:, ptt_idx_xy) * (XYmat_ortho(:, ptt_idx_xy)' * p);
p_orth = p - proj_on_ptt;

if norm(p_orth) > 1e-12
    p_orth = p_orth / norm(p_orth);
    removal_basis_xy = [XYmat_ortho(:, ptt_idx_xy), p_orth]; % 用于去除 PTT+power
    coeffs_obs = removal_basis_xy \ Wz;
    resid_pttpower_removed_xy = Wz - removal_basis_xy * coeffs_obs;
    RMS_XY_PTTPower_Removed = rms(resid_pttpower_removed_xy);
    PV_XY_PTTPower_Removed = max(resid_pttpower_removed_xy) - min(resid_pttpower_removed_xy);
   % fprintf('（XY）去除 PTT+Power 后 RMS=%.2f nm, PV=%.2f nm\n', RMS_XY_PTTPower_Removed*1e6, PV_XY_PTTPower_Removed*1e6);
    % 对拟合的 xy_fit，计算其在 removal_basis_xy 上的投影（以便画出去除后的拟合面）
    coeffs_fit = removal_basis_xy \ xy_fit;
    PTTPower_term_xy_fit = removal_basis_xy * coeffs_fit; % same unit
else
    resid_pttpower_removed_xy = Wz;
    PTTPower_term_xy_fit = zeros(size(Wz));
    warning('构��的 XY power 向量在与 PTT 正交后接近 0，跳过 XY PTT+Power 去除（拟合面）处理。');
end

% Zern_Fitting_Results=[RMS;length(zernike_coef);rms(zernike_resid) * 1e6;RMS_Zern_PTT_Removed*1e6;PV_Zern_PTT_Removed*1e6;RMS_Zern_PTTPower_Removed*1e6;PV_Zern_PTTPower_Removed*1e6];
% 
% Zern_Fitting_Iterm={'FEA输入表面误差RMS (nm)';'拟合项数';'拟合残差 (nm)';'去PPT项后RMS (nm)';'去PPT项后PV (nm)';'去PPTP项后RMS (nm)';'去PPTP项后PV (nm)'};
% Zern_Fitting_Matrix=table(Zern_Fitting_Iterm,Zern_Fitting_Results)
% 
% XY_Fitting_Results=[RMS;length(xy_coef);rms(xy_resid) * 1e6;RMS_XY_PTT_Removed*1e6; PV_XY_PTT_Removed*1e6;RMS_XY_PTTPower_Removed*1e6; PV_XY_PTTPower_Removed*1e6];
% 
% XY_Fitting_Iterm={'FEA输入表面误差RMS (nm)';'拟合项数';'拟合残差 (nm)';'去PPT项后RMS (nm)';'去PPT项后PV (nm)';'去PPTP项后RMS (nm)';'去PPTP项后PV (nm)'};
% XY_Fitting_Matrix=table(Fitting_Iterm,Fitting_Results)

Zern_Fitting_Results=[RMS;length(zernike_coef);rms(zernike_resid) * 1e6;RMS_Zern_PTT_Removed*1e6;PV_Zern_PTT_Removed*1e6;RMS_Zern_PTTPower_Removed*1e6;PV_Zern_PTTPower_Removed*1e6];

XY_Fitting_Results=[RMS;length(xy_coef);rms(xy_resid) * 1e6;RMS_XY_PTT_Removed*1e6; PV_XY_PTT_Removed*1e6;RMS_XY_PTTPower_Removed*1e6; PV_XY_PTTPower_Removed*1e6];

Fitting_Iterm={'FEA输入表面误差RMS (nm)';'拟合项数';'拟合残差 (nm)';'去PPT项后RMS (nm)';'去PPT项后PV (nm)';'去PPTP项后RMS (nm)';'去PPTP项后PV (nm)'};

Fitting_Matrix=table(Fitting_Iterm,Zern_Fitting_Results,XY_Fitting_Results);

Fitting_Matrix.Fitting_Iterm=char(Fitting_Iterm);
fprintf('\n\n===光学表面误差拟合===\n\n');
disp(Fitting_Matrix)

% XY_Fitting_Iterm={'FEA输入表面误差RMS (nm)';'拟合项数';'拟合残差 (nm)';'去PPT项后RMS (nm)';'去PPT项后PV (nm)';'去PPTP项后RMS (nm)';'去PPTP项后PV (nm)'};
% XY_Fitting_Matrix=table(Fitting_Iterm,Fitting_Results)

%%
%% ---------- 计算 slope RMS（grid 与 lsplane）并集成到输出 ----------
% 配置（可调整）
slope_opts_grid.method = 'grid';
slope_opts_grid.grid_nx = 200;        % 网格分辨率（grid 方法）
slope_opts_grid.remove_tilt = true;   % 是否去除最佳拟合平面（通常要去）
slope_opts_grid.interp_method = 'natural';

slope_opts_ls.method = 'lsplane';
slope_opts_ls.k_neighbors = 12;       % 局部邻域大小（lsplane 方法）
slope_opts_ls.remove_tilt = true;

% 输入：X,Y 单位 mm（与你脚本一致）；Wz 单位 m（脚本中 Wz = UZ，注：为 m）
% grid 方法返回基于插值网格的斜率 RMS（rad 单位）
[sx_rms_g, sy_rms_g, s_mag_rms_g, slopes_grid, info_grid] = compute_slope_rms(X, Y, W(:,3), slope_opts_grid);

% lsplane 方法在原始点上为每点拟合局部平面，返回按点的 sx, sy（rad 单位）及 RMS
[sx_rms_l, sy_rms_l, s_mag_rms_l, slopes_ls, info_ls] = compute_slope_rms(X, Y, W(:,3), slope_opts_ls);

% 将 slope RMS 转换为更常用单位（urad）
slopeRMS_grid_urad = s_mag_rms_g * 1e6;
slopeRMS_lsplane_urad = s_mag_rms_l * 1e6;

fprintf('Slope RMS (grid)    : %.3f urad\n', slopeRMS_grid_urad);
fprintf('Slope RMS (lsplane) : %.3f urad\n', slopeRMS_lsplane_urad);

% 将 slope RMS 写入 plot_pv_rms_summary.csv：在原有 summary 表中增加两列（每行重复，便于合并）
% 如果 summary_table 尚未构造（你的脚本在循环后才构造），这里我们把新列准备好并在写表时合并。
% 假设你已有 summary_table 生成代码如下（替换/合并时请确保列数一致）：
% summary_table = table(summary_idx, summary_title', summary_pv, summary_rms, summary_cmin, summary_cmax, ...
%    'VariableNames', {'Index','Title','PV_nm','RMS_nm','ColorMin_nm','ColorMax_nm'});

% 如果 summary_table 已在脚本后面创建，那么在写表之前执行下面两行将新列附加：
% 下面示例演示如何把新列填充为重复值，然后写出 CSV（在创建 summary_table 后运行）
% 例如（如果你之后有 writetable(summary_table, 'plot_pv_rms_summary.csv'); 这一句）：
% summary_table.SlopeRMS_grid_urad = repmat(slopeRMS_grid_urad, height(summary_table), 1);
% summary_table.SlopeRMS_lsplane_urad = repmat(slopeRMS_lsplane_urad, height(summary_table), 1);

% 如果你的脚本在此处已创建 summary_table（如示例脚本），下面将直接扩展并覆盖写文件：
if exist('summary_table','var')
    summary_table.SlopeRMS_grid_urad = repmat(slopeRMS_grid_urad, height(summary_table), 1);
    summary_table.SlopeRMS_lsplane_urad = repmat(slopeRMS_lsplane_urad, height(summary_table), 1);
    writetable(summary_table, 'plot_pv_rms_summary.csv');
    fprintf('已把 PV/RMS summary 与 SlopeRMS 列写入 plot_pv_rms_summary.csv\n');
else
    % 如果 summary_table 尚未创建，保存单独的 slope summary 文件
    slope_summary = table(slopeRMS_grid_urad, slopeRMS_lsplane_urad, ...
        'VariableNames', {'SlopeRMS_grid_urad','SlopeRMS_lsplane_urad'});
    writetable(slope_summary, 'slope_rms_summary.csv');
    fprintf('summary_table 尚未在工作区创建，已把 slope RMS 写入 slope_rms_summary.csv\n');
end

% 在 mirror_fit_result_sector.csv 中添加每节点的局部斜率（使用 lsplane 得到的点级估计）
% slopes_ls.sx / sy 存放每个原始点的 slope（rad），与原始节点顺序一致
if isfield(slopes_ls,'sx') && numel(slopes_ls.sx) == numel(nodeID)
    Sx_urad = slopes_ls.sx(:) * 1e6;
    Sy_urad = slopes_ls.sy(:) * 1e6;
    Smag_urad = sqrt(Sx_urad.^2 + Sy_urad.^2);
    % 如果 res_t 尚未生成，按你脚本原本的变量构造：
    if exist('res_t','var')
        res_t.SlopeX_urad = Sx_urad;
        res_t.SlopeY_urad = Sy_urad;
        res_t.SlopeMag_urad = Smag_urad;
        writetable(res_t, 'mirror_fit_result_sector.csv');
        fprintf('已把每节点的 SlopeX/Y/Mag 写入 mirror_fit_result_sector.csv\n');
    else
        % 重新构造类似原脚本的 res_t 并写出
        res_t = table(nodeID, X, Y, Wz*1e6, zernike_fit*1e6, zernike_resid*1e6, xy_fit*1e6, xy_resid*1e6, ...
            Sx_urad, Sy_urad, Smag_urad, ...
            'VariableNames', {'NodeID','X','Y','Zerr_nm','ZernikeFit_nm','ZernikeResid_nm','XYFit_nm','XYResid_nm', ...
            'SlopeX_urad','SlopeY_urad','SlopeMag_urad'});
        writetable(res_t, 'mirror_fit_result_sector.csv');
        fprintf('已创建并写入 mirror_fit_result_sector.csv（含 slope 列）\n');
    end
else
    warning('未能获得逐点斜率，跳过向 mirror_fit_result_sector.csv 写入 slope 列。');
end



%% 插值到规则网格（用于 surf）
nx = 200; ny = 200;
xgv = linspace(min(X), max(X), nx);
ygv = linspace(min(Y), max(Y), ny);
[Xg, Yg] = meshgrid(xgv, ygv);

F_orig = scatteredInterpolant(X, Y, Wz, 'natural', 'none');
F_zern = scatteredInterpolant(X, Y, zernike_fit, 'natural', 'none');
F_xy   = scatteredInterpolant(X, Y, xy_fit, 'natural', 'none');
F_zern_ptt = scatteredInterpolant(X, Y, zernike_fit - PTT_term_zern, 'natural', 'none');
F_xy_ptt   = scatteredInterpolant(X, Y, xy_fit - PTT_term_xy, 'natural', 'none');
F_zern_pttpower = scatteredInterpolant(X, Y, zernike_fit - PTTPower_term_zern, 'natural', 'none');
F_xy_pttpower   = scatteredInterpolant(X, Y, xy_fit - PTTPower_term_xy_fit, 'natural', 'none');
F_zern_resid = scatteredInterpolant(X, Y, zernike_resid, 'natural', 'none');
F_xy_resid   = scatteredInterpolant(X, Y, xy_resid, 'natural', 'none');

% 将插值结果转换为 nm（用于 color）
Z_orig_grid = F_orig(Xg, Yg) * 1e6;          % nm
Z_zern_grid = F_zern(Xg, Yg) * 1e6;          % nm
Z_xy_grid   = F_xy(Xg, Yg) * 1e6;            % nm
Z_zern_ptt_grid = F_zern_ptt(Xg, Yg) * 1e6;  % nm
Z_xy_ptt_grid   = F_xy_ptt(Xg, Yg) * 1e6;    % nm
Z_zern_pttpower_grid = F_zern_pttpower(Xg, Yg) * 1e6; % nm
Z_xy_pttpower_grid   = F_xy_pttpower(Xg, Yg) * 1e6;   % nm
Z_zern_resid_grid = F_zern_resid(Xg, Yg) * 1e6; % nm
Z_xy_resid_grid   = F_xy_resid(Xg, Yg) * 1e6;   % nm

% 计算统一主表面色阶
all_vals = [Z_orig_grid(:); Z_zern_grid(:); Z_xy_grid(:); Z_zern_ptt_grid(:); Z_xy_ptt_grid(:); Z_zern_pttpower_grid(:); Z_xy_pttpower_grid(:)];
all_vals = all_vals(~isnan(all_vals));
if ~isempty(all_vals)
    surf_cmin = prctile(all_vals, 1);
    surf_cmax = prctile(all_vals, 99);
else
    surf_cmin = [];
    surf_cmax = [];
end
% 残差色阶
res_vals = [Z_zern_resid_grid(:); Z_xy_resid_grid(:)];
res_vals = res_vals(~isnan(res_vals));
if ~isempty(res_vals)
    resid_cmin = prctile(res_vals, 1);
    resid_cmax = prctile(res_vals, 99);
else
    resid_cmin = [];
    resid_cmax = [];
end

% --- 自动计算几何高度压缩因子（仅用于视觉，不影响 color 映射） ---
% Xg,Yg 单位为原始单位（mm）; Z_*_grid 单位为 nm => 转回 mm 用于比例估计
Z_sample_mm = Z_orig_grid / 1e6; % mm
if isempty(Z_sample_mm(:)) || all(isnan(Z_sample_mm(:)))
    Z_range_mm = 0;
else
    Z_range_mm = max(Z_sample_mm(~isnan(Z_sample_mm))) - min(Z_sample_mm(~isnan(Z_sample_mm)));
end
XY_range_mm = max([max(Xg(:))-min(Xg(:)), max(Yg(:))-min(Yg(:))]);
if Z_range_mm == 0
    zscale = 1e-3;
else
    target_frac = 0.08; % 几何高度占 XY 范围的比例，可调整
    zscale = (XY_range_mm * target_frac) / Z_range_mm;
end
zscale = min(max(zscale, 1e-6), 1e3);
%%
% fig = figure('Name','扇形子镜表面拟合（auto-limits surf）','Units','normalized','Position',[0.03 0.03 0.94 0.92]);
% 
% % 1 原始
% ax1 = subplot(3,3,1, 'Parent', fig);
% surf_auto(ax1, Xg, Yg, Z_orig_grid, zscale, surf_cmin, surf_cmax);
% title('1. 原始表面误差');
% 
% % 2 Zernike 拟合
% ax2 = subplot(3,3,2, 'Parent', fig);
% surf_auto(ax2, Xg, Yg, Z_zern_grid, zscale, surf_cmin, surf_cmax);
% title('2. Zernike 拟合表面');
% 
% % 3 XY 拟合
% ax3 = subplot(3,3,3, 'Parent', fig);
% surf_auto(ax3, Xg, Yg, Z_xy_grid, zscale, surf_cmin, surf_cmax);
% title('3. XY 多项式 拟合表面');
% 
% % 4 Zernike 去 PTT
% ax4 = subplot(3,3,4, 'Parent', fig);
% surf_auto(ax4, Xg, Yg, Z_zern_ptt_grid, zscale, surf_cmin, surf_cmax);
% title('4. Zernike 去除 PTT');
% 
% % 5 XY 去 PTT
% ax5 = subplot(3,3,5, 'Parent', fig);
% surf_auto(ax5, Xg, Yg, Z_xy_ptt_grid, zscale,surf_cmin, surf_cmax);
% title('5. XY 去除 PTT');
% 
% % 6 Zernike 去 PTT+Power
% ax6 = subplot(3,3,6, 'Parent', fig);
% surf_auto(ax6, Xg, Yg, Z_zern_pttpower_grid, zscale, surf_cmin, surf_cmax);
% title('6. Zernike 去除 PTT+Power');
% 
% % 7 XY 去 PTT+Power
% ax7 = subplot(3,3,7, 'Parent', fig);
% surf_auto(ax7, Xg, Yg, Z_xy_pttpower_grid, zscale, surf_cmin, surf_cmax);
% title('7. XY 去除 PTT+Power');
% 
% % 8 Zernike 残差（独立色阶）
% ax8 = subplot(3,3,8, 'Parent', fig);
% surf_auto(ax8, Xg, Yg, Z_zern_resid_grid, zscale, resid_cmin, resid_cmax);
% title('8. Zernike 拟合 残差 (nm)');
% 
% % 9 XY 残差（独立色阶）
% ax9 = subplot(3,3,9, 'Parent', fig);
% surf_auto(ax9, Xg, Yg, Z_xy_resid_grid, zscale, resid_cmin, resid_cmax);
% title('9. XY 拟合 残差 (nm)');
% 
% % Improve appearance
% colormap(fig, parula);
% axs = [ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9];
% for k = 1:numel(axs)
%     axes(axs(k));
%     camlight('headlight');
%     material('dull');
% end
% set(findall(fig,'-property','FontSize'),'FontSize',10);
% 只影响surf可视化，自动每幅图都自适应Z轴和色条范围



% figure('Name','各曲面及残差（每幅图标尺自适应）','Units','normalized','Position',[0.05 0.05 0.8 0.8]);
% grids = {Z_orig_grid, Z_zern_grid, Z_xy_grid, Z_zern_ptt_grid, Z_xy_ptt_grid, ...
%          Z_zern_pttpower_grid, Z_xy_pttpower_grid, Z_zern_resid_grid, Z_xy_resid_grid};
% titles = {'1. 原始表面误差', ...
%           '2. Zernike 拟合表面', ...
%           '3. XY 多项式 拟合表面', ...
%           '4. Zernike 去除 PTT', ...
%           '5. XY 去除 PTT', ...
%           '6. Zernike 去除 PTT+Power', ...
%           '7. XY 去除 PTT+Power', ...
%           '8. Zernike 拟合 残差 (nm)', ...
%           '9. XY 拟合 残差 (nm)'};
% for k = 1:9
%     ax = subplot(3,3,k);
%     Zgrid = grids{k};
%     surf(Xg, Yg, Zgrid, 'EdgeColor','none');
%     xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Surface (nm)');
%     view(2); axis equal; axis tight; shading interp; colorbar;
%     colormap(jet);
%     title(titles{k},'Interpreter','none');
%     zvals = Zgrid(~isnan(Zgrid));
%     if ~isempty(zvals)
%         zmin = min(zvals);
%         zmax = max(zvals);
% %         if zmin==zmax % 极小变形量时微调
% %             delta = max(abs(zmin)*0.1, 1e-8);
% %             zlim([zmin-delta zmax+delta]);
% %             caxis([zmin-delta zmax+delta]);
% %         else
% %             zlim([zmin zmax]);
%             caxis([zmin zmax]);
% %         end
%     end
%     camlight('headlight');
%     material('dull');
% end

% 可视化替换块：使用 UZ 作为数据源，按百分位剪切 caxis，
% 在右上角显示 PV / RMS（nm），并导出 summary CSV。
% 将此段替换原脚本中插值并绘图的部分（仅替换可视化，不更改拟合计算）

% 预期：脚本在此之前已计算好网格 Xg,Yg（单位 mm）
% 并计算好以下 Z_*_grid（单位 nm）：
% Z_orig_grid, Z_zern_grid, Z_xy_grid, Z_zern_ptt_grid, Z_xy_ptt_grid,
% Z_zern_pttpower_grid, Z_xy_pttpower_grid, Z_zern_resid_grid, Z_xy_resid_grid

% 配置（可按需调整）
p_low = 1;   % lower percentile for clipping (percent)
p_high = 99; % upper percentile for clipping (percent)
colormap_choice = @jet; % 颜色图（最大值为红色）

% 以 cell 组织要绘制的数据与标题（按子图顺序）
% grids = {Z_orig_grid, Z_zern_grid, Z_xy_grid, Z_zern_ptt_grid, Z_xy_ptt_grid, ...
%          Z_zern_pttpower_grid, Z_xy_pttpower_grid, Z_zern_resid_grid, Z_xy_resid_grid};
% titles = {'1. 原始表面误差', ...
%           '2. Zernike拟合表面', ...
%           '3. XY 多项式拟合表面', ...
%           '4. Zernike去除PTT(Piston&Tip&Tilt)', ...
%           '5. XY去除去除PTT', ...
%           '6. Zernike去除PTTP(Piston&Tip&Tilt&Power)', ...
%           '7. XY去除PTTP', ...
%           '8. Zernike拟合残差', ...
%           '9. XY拟合残差'};
grids = {Z_orig_grid, Z_zern_grid, Z_xy_grid, Z_zern_ptt_grid, Z_zern_pttpower_grid, ...
         Z_zern_resid_grid,Z_xy_ptt_grid, Z_xy_pttpower_grid, Z_xy_resid_grid};
titles = {'1. 原始表面误差', ...
          '2. Zernike拟合表面', ...
          '3. XY多项式拟合表面', ...
          '4. Zernike去除PTT(Piston&Tip&Tilt)', ...
          '5. Zernike去除PTTP(Piston&Tip&Tilt&Power)', ...
          '6. Zernike拟合残差', ...
          '7. XY去除PTT', ...
          '8. XY去除PTTP', ...
          '9. XY拟合残差'};
% 输出 summary 初始化
nplots = numel(grids);
summary_idx = (1:nplots).';
summary_title = titles;
summary_pv = nan(nplots,1);
summary_rms = nan(nplots,1);
summary_cmin = nan(nplots,1);
summary_cmax = nan(nplots,1);

% 创建图像窗口
fig = figure('Name','扇形子镜表面拟合（带 PV/RMS 与百分位剪切）','Units','normalized','Position',[0.03 0.03 0.8 0.8]);
colormap(colormap_choice());

% 遍历绘图
% for k = 1:nplots
%     Zgrid = grids{k}; % 单位：nm
%     ax = subplot(3,3,k, 'Parent', fig);
%     
%     % 统计有效值（排除 NaN/Inf）
%     valid = Zgrid(~isnan(Zgrid) & ~isinf(Zgrid));
%     if isempty(valid)
%         title(titles{k}, 'Interpreter','none');
%         warning('第 %d 幅图 "%s" 无有效数据，跳过绘制。', k, titles{k});
%         continue;
%     end
%     
%     % 计算 PV 与 RMS（基于有效数据，单位 nm）
%     pv_val = max(valid) - min(valid);
%     rms_val = rms(valid);
%     summary_pv(k) = pv_val;
%     summary_rms(k) = rms_val;
%     
%     % 计算 percentile 剪切范围（用于 color 映射）
%     cmin = prctile(valid, p_low);
%     cmax = prctile(valid, p_high);
%     % 如果 p_low/p_high 恰好使 cmin==cmax（极少见），则退回 min/max
%     if cmin == cmax
%         cmin = min(valid);
%         cmax = max(valid);
%     end
%     summary_cmin(k) = cmin;
%     summary_cmax(k) = cmax;
%     
%     % 为避免 colorbar 被少数极端值影响，采用 cmin/cmax 映射色阶；几何 zlim 使用真实 min/max
%     zmin_all = min(valid); % 真正的 min (nm)
%     zmax_all = max(valid); % 真正的 max (nm)
%     
%     % 为几何显示计算缩放（仍保留，使视觉不被 X/Y 数值范围和 Z 的数量级差异干扰）
%     % 将 nm -> mm 用作几何高度基础
%     Z_mm = Zgrid / 1e6; % mm
%     % 计算 XY 范围（mm）
%     XY_range = max(max(Xg(:)) - min(Xg(:)), max(Yg(:)) - min(Yg(:)));
%     % Z 范围（mm）
%     if zmax_all == zmin_all
%         zrange_mm = max(abs(zmin_all),1) / 1e6;
%     else
%         zrange_mm = (zmax_all - zmin_all) / 1e6;
%     end
%     % target fraction of XY range to occupy by Z geometry
%     target_frac = 0.08;
%     if zrange_mm == 0
%         zscale = 1e-3;
%     else
%         zscale = (XY_range * target_frac) / zrange_mm;
%         zscale = min(max(zscale, 1e-6), 1e3);
%     end
%     
%     % 生成几何高度（mm），同时颜色映射使用原始 nm 值
%     Zgeom = Z_mm * zscale;
%     hsurf = surf(Xg, Yg, Zgeom, 'EdgeColor','none');
%     set(hsurf, 'CData', Zgrid); % 使用 nm 数据作为 CData -> colorbar 显示 nm
%     shading interp;
%     view(2);
%     axis equal; axis tight;
%     xlabel('X (mm)'); ylabel('Y (mm)');
%     zlabel('Surface (nm)');
%     title(titles{k}, 'Interpreter','none');
%     
%     % 设置 colorbar（使用 percentile 剪切）
%     cb = colorbar;
%     caxis([cmin cmax]);  % color 映射使用百分位剪切范围（nm）
%     ylabel(cb, 'Surface (nm)');
%     
%     % 设置几何 zlim 映射为数据真实 min/max（转换为几何高度后）
%     vis_zmin = (zmin_all / 1e6) * zscale;
%     vis_zmax = (zmax_all / 1e6) * zscale;
%     if vis_zmin == vis_zmax
%         % 微调避免单点区间
%         delta = max(abs(vis_zmax)*0.01, 1e-6);
%         zlim([vis_zmin-delta, vis_zmax+delta]);
%     else
%         zlim([vis_zmin, vis_zmax]);
%     end
%     
%     % 将几何的 Z 轴刻度标签映射为 nm 显示（更直观）
%     zticks_vis = get(gca, 'ZTick');
%     zticklabels_nm = arrayfun(@(v) sprintf('%.0f', (v / zscale) * 1e6), zticks_vis, 'UniformOutput', false);
%     set(gca, 'ZTickLabel', zticklabels_nm);
%     
%     % 将 PV/RMS 添加到右上角（normalized 轴坐标）
%     txt = sprintf('PV=%.2f nm\nRMS=%.2f nm', pv_val, rms_val);
%     % 使用 annotation 风格也可，但这里用轴内 text（normalized）
%     text(0.98, 0.98, txt, ...
%         'Units', 'normalized', ...
%         'HorizontalAlignment', 'right', ...
%         'VerticalAlignment', 'top', ...
%         'FontSize', 10, ...
%         'FontWeight', 'bold', ...
%         'Color',[0,0,0],...
%         'BackgroundColor', 'none', ...
%         'EdgeColor', 'none', ...
%         'Interpreter', 'none');
%     
%     % 光照和材质
%     camlight('headlight');
%     material('dull');
% end
% 
% % 导出 PV/RMS summary 到 CSV
% summary_table = table(summary_idx, summary_title', summary_pv, summary_rms, summary_cmin, summary_cmax, ...
%     'VariableNames', {'Index','Title','PV_nm','RMS_nm','ColorMin_nm','ColorMax_nm'});
% writetable(summary_table, 'plot_pv_rms_summary.csv');
% fprintf('PV/RMS summary 已保存到 plot_pv_rms_summary.csv\n');
% 
% % 美化字体
% set(findall(fig,'-property','FontSize'),'FontSize',10);
% 
% %% 结果导出（保持原样）
% res_t = table(nodeID, X, Y, Wz*1e6, zernike_fit*1e6, zernike_resid*1e6, xy_fit*1e6, xy_resid*1e6, ...
%     'VariableNames', {'NodeID','X','Y','Zerr_nm','ZernikeFit_nm','ZernikeResid_nm','XYFit_nm','XYResid_nm'});
% writetable(res_t, 'mirror_fit_result_sector.csv');
% writematrix(zernike_coef, 'zernike_coefs_sector.csv');
% writematrix(xy_coef, 'xy_poly_coefs_sector.csv');
% 
% fprintf('绘图完成：每幅图的坐标范围已依据实际数据自动设置（X/Y 单位保留原始，Z 以 nm 显示）。\n');
% 
% 
% % Helper: 绘制并自动设置每图 x/y/z limits 与 ztick labels 映射到 nm
% 
% 
% 
% function surf_auto(ax, Xg, Yg, Z_nm, zscale, cmin, cmax)
%     axes(ax);
%     % 几何高度（mm）
%     Z_vis = (Z_nm / 1e6) * zscale; % mm
%     hsurf = surf(Xg, Yg, Z_vis, 'EdgeColor','none');
%     set(hsurf, 'CData', Z_nm); % 颜色代表 nm
%     shading interp; view(2); axis equal; axis tight;
%     xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Surface (nm)');
%     % 设置 xy limits 为网格范围
%     xlim([min(Xg(:)), max(Xg(:))]); ylim([min(Yg(:)), max(Yg(:))]);
%     % 计算并设置 z limits（基于真实 Z_nm）
%     Zvals = Z_nm(:); Zvals = Zvals(~isnan(Zvals));
%     if ~isempty(Zvals)
%         zmin_nm = min(Zvals); zmax_nm = max(Zvals);
%         zmin_vis = (zmin_nm / 1e6) * zscale;
%         zmax_vis = (zmax_nm / 1e6) * zscale;
%         if zmin_vis == zmax_vis
%             % 给出一个小的宽度
%             delta = max(abs(zmin_vis)*0.01, 1e-6);
%             zlim([zmin_vis-delta, zmax_vis+delta]);
%         else
%             zlim([zmin_vis, zmax_vis]);
%         end
%         % 设置 z ticks (5 ticks) 并把标签映为 nm
%         zticks_vis = linspace(zmin_vis, zmax_vis, 5);
%         set(gca, 'ZTick', zticks_vis);
%         zticklabels_nm = arrayfun(@(v) sprintf('%.0f', (v / zscale) * 1e6), zticks_vis, 'UniformOutput', false);
%         set(gca, 'ZTickLabel', zticklabels_nm);
%     end
%     cb = colorbar;
%     ylabel(cb, 'Surface (nm)');
%     if ~isempty(cmin) && ~isempty(cmax)
%         caxis([cmin cmax]);
%     end
%     grid on;
% end
% 
% % 绘图（3x3）
% 替换：for k = 1:nplots ... end 循环体（绘图与 PV/RMS 标题）
for k = 1:nplots
    Zgrid = grids{k}; % 单位：nm
    ax = subplot(3,3,k, 'Parent', fig);
    
    % 统计有效值（排除 NaN/Inf）
    valid = Zgrid(~isnan(Zgrid) & ~isinf(Zgrid));
    if isempty(valid)
        % 如果无有效数据，仍显示标题（仅标题）
        title({titles{k}, 'PV: -- nm    RMS: -- nm'}, 'Interpreter', 'none');
        warning('第 %d 幅图 "%s" 无有效数据，跳过绘制。', k, titles{k});
        continue;
    end
    
    % 计算 PV 与 RMS（基于有效数据，单位 nm）
    pv_val = max(valid) - min(valid);
    rms_val = rms(valid);
    summary_pv(k) = pv_val;
    summary_rms(k) = rms_val;
    
    % 计算 percentile 剪切范围（用于 color 映射）
    cmin = prctile(valid, p_low);
    cmax = prctile(valid, p_high);
    % 如果 p_low/p_high 恰好使 cmin==cmax（极少见），则退回 min/max
    if cmin == cmax
        cmin = min(valid);
        cmax = max(valid);
    end
    summary_cmin(k) = cmin;
    summary_cmax(k) = cmax;
    
    % 为避免 colorbar 被少数极端值影响，采用 cmin/cmax 映射色阶；几何 zlim 使用真实 min/max
    zmin_all = min(valid); % 真正的 min (nm)
    zmax_all = max(valid); % 真正的 max (nm)
    
    % 为几何显示计算缩放（仍保留，使视觉不被 X/Y 数值范围和 Z 的数量级差异干扰）
    % 将 nm -> mm 用作几何高度基础
    Z_mm = Zgrid / 1e6; % mm
    % 计算 XY 范围（mm）
    XY_range = max(max(Xg(:)) - min(Xg(:)), max(Yg(:)) - min(Yg(:)));
    % Z 范围（mm）
    if zmax_all == zmin_all
        zrange_mm = max(abs(zmin_all),1) / 1e6;
    else
        zrange_mm = (zmax_all - zmin_all) / 1e6;
    end
    % target fraction of XY range to occupy by Z geometry
    target_frac = 0.08;
    if zrange_mm == 0
        zscale = 1e-3;
    else
        zscale = (XY_range * target_frac) / zrange_mm;
        zscale = min(max(zscale, 1e-6), 1e3);
    end
    
    % 生成几何高度（mm），同时颜色映射使用原始 nm 值
    Zgeom = Z_mm * zscale;
    hsurf = surf(Xg, Yg, Zgeom, 'EdgeColor','none');
    set(hsurf, 'CData', Zgrid); % 使用 nm 数据作为 CData -> colorbar 显示 nm
    shading interp;
    view(2);
    axis equal; axis tight;
    xlabel('X (mm)'); ylabel('Y (mm)');
    zlabel('Surface (nm)');
    
    % 使用两行标题：第一行为原始标题，第二行为 PV / RMS（位于标题下方，图像上方）
    % 使用 cell array 可以创建多行标题。'Interpreter','none' 保持原始文本不被解释。
    title({titles{k}, sprintf('PV = %.2f nm    RMS = %.2f nm', pv_val, rms_val)}, 'Interpreter', 'none');
    
    % 设置 colorbar（使用 percentile 剪切）
    cb = colorbar;
    caxis([cmin cmax]);  % color 映射使用百分位剪切范围（nm）
    ylabel(cb, 'Surface (nm)');
    
    % 设置几何 zlim 映射为数据真实 min/max（转换为几何高度后）
    vis_zmin = (zmin_all / 1e6) * zscale;
    vis_zmax = (zmax_all / 1e6) * zscale;
    if vis_zmin == vis_zmax
        % 微调避免单点区间
        delta = max(abs(vis_zmax)*0.01, 1e-6);
        zlim([vis_zmin-delta, vis_zmax+delta]);
    else
        zlim([vis_zmin, vis_zmax]);
    end
    
    % 将几何的 Z 轴刻度标签映射为 nm 显示（更直观）
    zticks_vis = get(gca, 'ZTick');
    zticklabels_nm = arrayfun(@(v) sprintf('%.0f', (v / zscale) * 1e6), zticks_vis, 'UniformOutput', false);
    set(gca, 'ZTickLabel', zticklabels_nm);
    
    % 光照和材质
    camlight('headlight');
    material('dull');
end
%%
% 替换：surf_auto 辅助函数（如果你仍然使用 surf_auto，应同步更新）
function surf_auto(ax, Xg, Yg, Z_nm, zscale, cmin, cmax, pv_val, rms_val)
    % 如果调用时没有传 pv_val, rms_val，可用 placeholder 文本
    if nargin < 8
        pv_val = [];
        rms_val = [];
    end

    axes(ax);
    % 几何高度（mm）
    Z_vis = (Z_nm / 1e6) * zscale; % mm
    hsurf = surf(Xg, Yg, Z_vis, 'EdgeColor','none');
    set(hsurf, 'CData', Z_nm); % 颜色代表 nm
    shading interp; view(2); axis equal; axis tight;
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Surface (nm)');
    % 设置 xy limits 为网格范围
    xlim([min(Xg(:)), max(Xg(:))]); ylim([min(Yg(:)), max(Yg(:))]);
    % 计算并设置 z limits（基于真实 Z_nm）
    Zvals = Z_nm(:); Zvals = Zvals(~isnan(Zvals));
    if ~isempty(Zvals)
        zmin_nm = min(Zvals); zmax_nm = max(Zvals);
        zmin_vis = (zmin_nm / 1e6) * zscale;
        zmax_vis = (zmax_nm / 1e6) * zscale;
        if zmin_vis == zmax_vis
            % 给出一个小的宽度
            delta = max(abs(zmin_vis)*0.01, 1e-6);
            zlim([zmin_vis-delta, zmax_vis+delta]);
        else
            zlim([zmin_vis, zmax_vis]);
        end
        % 设置 z ticks (5 ticks) 并把标签映为 nm
        zticks_vis = linspace(zmin_vis, zmax_vis, 5);
        set(gca, 'ZTick', zticks_vis);
        zticklabels_nm = arrayfun(@(v) sprintf('%.0f', (v / zscale) * 1e6), zticks_vis, 'UniformOutput', false);
        set(gca, 'ZTickLabel', zticklabels_nm);
    end
    cb = colorbar;
    ylabel(cb, 'Surface (nm)');
    if ~isempty(cmin) && ~isempty(cmax)
        caxis([cmin cmax]);
    end
    % 标题：两行（若有 pv/rms 则显示，否则用占位）
    if ~isempty(pv_val) && ~isempty(rms_val)
        ttext = sprintf('PV = %.2f nm    RMS = %.2f nm', pv_val, rms_val);
    else
        ttext = 'PV: -- nm    RMS: -- nm';
    end
    % 注意：调用 surf_auto 时传入的标题应由调用者设置；这里仅展示 pv/rms 可选的两行标题
    % 如果希望同时显示主标题与 pv/rms，请在外部用 title({mainTitle, ttext}, 'Interpreter','none')
    grid on;
end
%% ========== Zernike基函数 ==========
function Zmat = generateZernikeBasis(rho, theta, order)
    N = length(rho);
    modes = [];
    for n = 0:order
        for m = -n:2:n
            modes = [modes; n, m];
        end
    end
    n_modes = size(modes, 1);
    Zmat = zeros(N, n_modes);
    for k = 1:n_modes
        n = modes(k,1); m = modes(k,2);
        Zvec = zeros(N,1);
        for i = 1:N
            if rho(i) <= 1
                Rnm = zernikeRadial(rho(i), n, abs(m));
                if m > 0
                    Zvec(i) = Rnm * cos(m * theta(i));
                elseif m < 0
                    Zvec(i) = Rnm * sin(-m * theta(i));
                else
                    Zvec(i) = Rnm;
                end
            else
                Zvec(i) = 0;
            end
        end
        fac = sqrt(sum(Zvec.^2));
        if fac > 1e-10
            Zvec = Zvec / fac;
        end
        Zmat(:,k) = Zvec;
    end
end

function Rnm = zernikeRadial(r, n, m)
    if mod(n - m, 2) ~= 0
        Rnm = 0; return;
    end
    Rnm = 0;
    upper = (n - m) / 2;
    for k = 0:upper
        coeff = ((-1)^k) * factorial(n - k) / (factorial(k) * factorial((n + m)/2 - k) * factorial((n - m)/2 - k));
        Rnm = Rnm + coeff * r^(n - 2*k);
    end
end

%% ========== XY多项式基函数 ==========
function [XYmat, terms] = generateXYPolynomialBasis(x, y, order)
    x = x(:); y = y(:);
    N = length(x); columns = {};
    terms = [];
    for i = 0:order
        for j = 0:order
            if i + j <= order
                columns{end+1} = x.^i .* y.^j;
                terms = [terms; i, j];
            end
        end
    end
    XYmat = [columns{:}];
    for k = 1:size(XYmat,2)
        normf = norm(XYmat(:,k));
        if normf > 1e-12
            XYmat(:,k) = XYmat(:,k) / normf;
        end
    end
end
%%
%% ====== 计算 slope 的函数（放在脚本末尾或独立文件 compute_slope_rms.m） ======
function [sx_rms, sy_rms, s_mag_rms, slopes, info] = compute_slope_rms(X, Y, Z_m, opts)
    % 输入：
    %   X, Y: 列向量，单位 mm（与你脚本一致）
    %   Z_m : 列向量，单位 m（Wz 在你的脚本中为 m）
    %   opts: 参数
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'method'), opts.method = 'grid'; end
    if ~isfield(opts,'grid_nx'), opts.grid_nx = 200; end
    if ~isfield(opts,'remove_tilt'), opts.remove_tilt = true; end
    if ~isfield(opts,'k_neighbors'), opts.k_neighbors = 12; end
    if ~isfield(opts,'interp_method'), opts.interp_method = 'natural'; end

    % 将 X,Y mm -> m
    Xv = X(:); Yv = Y(:); Zv = Z_m(:);
    X_m = Xv ;
    Y_m = Yv ;
    Zm = Zv; % 已是 m

    info = opts;

    % 去除最佳拟合平面（通常用于移除全局 tilt）
    if opts.remove_tilt
        A = [X_m, Y_m, ones(size(X_m))];
        coeff = A \ Zm;
        Zm = Zm - (A * coeff);
        info.removed_plane = coeff;
    else
        info.removed_plane = [0;0;0];
    end

    switch lower(opts.method)
        case 'grid'
            nx = opts.grid_nx;
            xgv = linspace(min(X_m), max(X_m), nx);
            ygv = linspace(min(Y_m), max(Y_m), nx);
            [Xg, Yg] = meshgrid(xgv, ygv);
            F = scatteredInterpolant(X_m, Y_m, Zm, opts.interp_method, 'none');
            Zg = F(Xg, Yg); % m
            valid_mask = ~isnan(Zg) & ~isinf(Zg);
            if nnz(valid_mask) == 0
                error('插值后无有效网格点，请改用 lsplane 方法或调整 grid_nx/插值方法。');
            end
            dx = xgv(2)-xgv(1);
            dy = ygv(2)-ygv(1);
            [dZdy, dZdx] = gradient(Zg, dy, dx); % 注意顺序：gradient(Z, dy, dx)
            sx_all = dZdx; sy_all = dZdy; % 单位：m/m = rad (小角近似)
            sx = sx_all(valid_mask); sy = sy_all(valid_mask);
            s_mag = sqrt(sx.^2 + sy.^2);
            sx_rms = sqrt(mean(sx.^2));
            sy_rms = sqrt(mean(sy.^2));
            s_mag_rms = sqrt(mean(s_mag.^2));
            slopes.sx = sx_all; slopes.sy = sy_all; slopes.mask = valid_mask;
            slopes.Xg = Xg; slopes.Yg = Yg; slopes.Zg = Zg;
            info.method = 'grid';
            info.grid_dx = dx; info.grid_dy = dy;
        case 'lsplane'
            pts = [X_m, Y_m];
            M = size(pts,1);
            K = min(opts.k_neighbors, M);
            % 使用 knnsearch 找邻居
            idxs = knnsearch(pts, pts, 'K', K);
            a_vec = zeros(M,1); b_vec = zeros(M,1);
            for i = 1:M
                nb = idxs(i,:);
                Xi = X_m(nb); Yi = Y_m(nb); Zi = Zm(nb);
                G = [Xi, Yi, ones(numel(Xi),1)];
                c = G \ Zi; % [a; b; c]
                a_vec(i) = c(1); b_vec(i) = c(2);
            end
            sx = a_vec; sy = b_vec; % 单位 rad
            s_mag = sqrt(sx.^2 + sy.^2);
            sx_rms = sqrt(mean(sx.^2));
            sy_rms = sqrt(mean(sy.^2));
            s_mag_rms = sqrt(mean(s_mag.^2));
            slopes.sx = sx; slopes.sy = sy; slopes.s_mag = s_mag;
            info.method = 'lsplane';
            info.k_neighbors = K;
        otherwise
            error('不支持的 method: %s', opts.method);
    end

    % 返回 info 中添加常用单位
    info.sx_rms_rad = sx_rms;
    info.sx_rms_urad = sx_rms * 1e6;
    info.sy_rms_rad = sy_rms;
    info.s_mag_rms_rad = s_mag_rms;
    info.s_mag_rms_urad = s_mag_rms * 1e6;
end