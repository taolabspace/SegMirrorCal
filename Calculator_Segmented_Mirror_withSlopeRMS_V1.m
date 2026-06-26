clear; clc; close all;

%% SigFit 对标版：扇形子镜光学表面误差计算
% 说明：
% 1) 输入 txt 文件第 2-7 列单位均为 m；脚本内部统一转换为 mm，输出 RMS/PV 为 nm。
% 2) 单独输出两套更接近 SigFit 的结果：
%    - Direct fit remove PTT：直接拟合下去除 piston/tip/tilt
%    - Best-fit plane remove PTT：最佳吻合面下去除 piston/tip/tilt
% 3) 保留 Zernike / XY 高阶拟合作为分析用途，但不再将其 full residual 作为 SigFit 对标主指标。

%% 0. 参数
input_file = '2.txt';
zernike_order = 10;
xy_order = 10;
remove_rigid_3d = false;   % 对标 SigFit 时建议关闭；true 会更强地吃掉低阶量

%% 1. 数据导入和单位转换（m -> mm）
data = readmatrix(input_file);

nodeID = data(:,1);
X  = data(:,2) * 1000;
Y  = data(:,3) * 1000;
Z  = data(:,4) * 1000;
UX = data(:,5) * 1000;
UY = data(:,6) * 1000;
UZ = data(:,7) * 1000;

orig_points = [X, Y, Z];
U = [UX, UY, UZ];
deformed_points = orig_points + U;

%% 2. 原始光学面局部坐标（用于 direct-fit）
orig_center = mean(orig_points, 1);
A0 = orig_points - orig_center;
[~, ~, V0] = svd(A0, 0);
ex0 = V0(:,1);
ey0 = V0(:,2);
n0  = V0(:,3);
n0  = n0 / norm(n0);

local_orig = A0 * [ex0 ey0 n0];
x0 = local_orig(:,1);
y0 = local_orig(:,2);

%% 3. Direct fit remove PTT（更接近“直接拟合下去除前三阶”）
% 沿原始光学面法向投影位移
W_direct = U * n0;   % mm
A_ptt_direct = [ones(size(x0)), x0, y0];
coef_direct_ptt = A_ptt_direct \ W_direct;
fit_direct_ptt = A_ptt_direct * coef_direct_ptt;
resid_direct_ptt = W_direct - fit_direct_ptt;

RMS_Direct_PTT_nm = rms(resid_direct_ptt) * 1e6;
PV_Direct_PTT_nm  = (max(resid_direct_ptt) - min(resid_direct_ptt)) * 1e6;

%% 4. 可选：三维刚体配准（默认关闭，仅供研究）
if remove_rigid_3d
    center_orig = mean(orig_points,1);
    center_def  = mean(deformed_points,1);
    P1 = orig_points - center_orig;
    P2 = deformed_points - center_def;
    H = P1' * P2;
    [U_svd, ~, V_svd] = svd(H);
    R = V_svd * U_svd';
    if det(R) < 0
        V_svd(:,3) = -V_svd(:,3);
        R = V_svd * U_svd';
    end
    T = center_def' - R * center_orig';
    P_corr = (R * orig_points' + T)';
    Residual = deformed_points - P_corr;
else
    R = eye(3);
    T = [0;0;0];
    P_corr = orig_points;
    Residual = U;
end

%% 5. Best-fit plane（更接近“最佳吻合面下剔除前三阶”）
% 用变形后点云求最佳吻合平面，再在该局部坐标系中去 PTT
center_def = mean(deformed_points, 1);
Ad = deformed_points - center_def;
[~, ~, Vd] = svd(Ad, 0);
exd = Vd(:,1);
eyd = Vd(:,2);
nd  = Vd(:,3);
nd  = nd / norm(nd);

local_def = Ad * [exd eyd nd];
xd = local_def(:,1);
yd = local_def(:,2);
wd = local_def(:,3);   % 相对最佳吻合面的 sag，单位 mm

A_ptt_bestfit = [ones(size(xd)), xd, yd];
coef_bestfit_ptt = A_ptt_bestfit \ wd;
fit_bestfit_ptt = A_ptt_bestfit * coef_bestfit_ptt;
resid_bestfit_ptt = wd - fit_bestfit_ptt;

RMS_BestFit_PTT_nm = rms(resid_bestfit_ptt) * 1e6;
PV_BestFit_PTT_nm  = (max(resid_bestfit_ptt) - min(resid_bestfit_ptt)) * 1e6;

%% 6. 统一高阶分析基准：使用 best-fit plane 局部坐标 + best-fit sag
% 注意：以下结果用于分析，不建议直接与 SigFit “去前三阶”列比较。
Xfit = xd;
Yfit = yd;
Wz   = wd;   % mm

PV = (max(Wz) - min(Wz)) * 1e6;
RMS = rms(Wz) * 1e6;

[roll, pitch, yaw] = dcm2angle(R, "XYZ");
RigidBody_Motions = [T(1)*1e3; T(2)*1e3; T(3)*1e3; roll*1e6; pitch*1e6; yaw*1e6];
RigidBody_Iterm = {'Tx (um)';'Ty (um)';'Tz (um)';'Rx (urad)';'Ry (urad)';'Rz (urad)'};
RigidBody_Matrix = table(RigidBody_Iterm, RigidBody_Motions);
RigidBody_Matrix.RigidBody_Iterm = char(RigidBody_Matrix.RigidBody_Iterm);
fprintf('\n===光学表面刚体运动（研究用）===\n\n');
disp(RigidBody_Matrix)

%% 7. Zernike 多项式拟合（分析用）
rho = hypot(Xfit, Yfit);
if max(rho) > 0
    rho = rho / max(rho);
end
theta = atan2(Yfit, Xfit);

Zmat = generateZernikeBasis(rho, theta, zernike_order);
[Zmat_ortho, ~] = qr(Zmat, 0);
zernike_coef = Zmat_ortho \ Wz;
zernike_fit = Zmat_ortho * zernike_coef;
zernike_resid = Wz - zernike_fit;

modes = [];
for n = 0:zernike_order
    for m = -n:2:n
        modes = [modes; n, m];
    end
end
ptt_idx = 1:min(3, size(Zmat_ortho,2));
defocus_idx = find(modes(:,1)==2 & modes(:,2)==0, 1);

PTT_term_zern = Zmat_ortho(:, ptt_idx) * zernike_coef(ptt_idx);
resid_ptt_removed_zern = Wz - PTT_term_zern;
RMS_Zern_PTT_Removed = rms(resid_ptt_removed_zern);
PV_Zern_PTT_Removed = max(resid_ptt_removed_zern) - min(resid_ptt_removed_zern);

if ~isempty(defocus_idx)
    idx_ptt_power = unique([ptt_idx(:); defocus_idx(:)]);
    PTTPower_term_zern = Zmat_ortho(:, idx_ptt_power) * zernike_coef(idx_ptt_power);
    resid_pttpower_removed_zern = Wz - PTTPower_term_zern;
    RMS_Zern_PTTPower_Removed = rms(resid_pttpower_removed_zern);
    PV_Zern_PTTPower_Removed = max(resid_pttpower_removed_zern) - min(resid_pttpower_removed_zern);
else
    PTTPower_term_zern = zeros(size(Wz));
    resid_pttpower_removed_zern = Wz;
    RMS_Zern_PTTPower_Removed = rms(resid_pttpower_removed_zern);
    PV_Zern_PTTPower_Removed = max(resid_pttpower_removed_zern) - min(resid_pttpower_removed_zern);
end

%% 8. XY 多项式拟合（分析用）
[XYmat, ~] = generateXYPolynomialBasis(Xfit, Yfit, xy_order);
[XYmat_ortho, ~] = qr(XYmat, 0);
xy_coef = XYmat_ortho \ Wz;
xy_fit = XYmat_ortho * xy_coef;
xy_resid = Wz - xy_fit;

ptt_idx_xy = 1:min(3, size(XYmat_ortho,2));
PTT_term_xy = XYmat_ortho(:, ptt_idx_xy) * xy_coef(ptt_idx_xy);
resid_ptt_removed_xy = Wz - PTT_term_xy;
RMS_XY_PTT_Removed = rms(resid_ptt_removed_xy);
PV_XY_PTT_Removed = max(resid_ptt_removed_xy) - min(resid_ptt_removed_xy);

power_vec = Xfit.^2 + Yfit.^2;
p = power_vec(:);
proj_on_ptt = XYmat_ortho(:, ptt_idx_xy) * (XYmat_ortho(:, ptt_idx_xy)' * p);
p_orth = p - proj_on_ptt;
if norm(p_orth) > 1e-12
    p_orth = p_orth / norm(p_orth);
    removal_basis_xy = [XYmat_ortho(:, ptt_idx_xy), p_orth];
    coeffs_obs = removal_basis_xy \ Wz;
    resid_pttpower_removed_xy = Wz - removal_basis_xy * coeffs_obs;
    RMS_XY_PTTPower_Removed = rms(resid_pttpower_removed_xy);
    PV_XY_PTTPower_Removed = max(resid_pttpower_removed_xy) - min(resid_pttpower_removed_xy);
    coeffs_fit = removal_basis_xy \ xy_fit;
    PTTPower_term_xy_fit = removal_basis_xy * coeffs_fit;
else
    resid_pttpower_removed_xy = Wz;
    RMS_XY_PTTPower_Removed = rms(resid_pttpower_removed_xy);
    PV_XY_PTTPower_Removed = max(resid_pttpower_removed_xy) - min(resid_pttpower_removed_xy);
    PTTPower_term_xy_fit = zeros(size(Wz));
    warning('构造的 XY power 向量在与 PTT 正交后接近 0，跳过 XY PTT+Power 去除。');
end

%% 9. 对标结果输出
fprintf('\n=== SigFit 对标主结果 ===\n\n');
Compare_Item = {
    'Direct fit remove PTT RMS (nm)';
    'Direct fit remove PTT PV (nm)';
    'Best-fit plane remove PTT RMS (nm)';
    'Best-fit plane remove PTT PV (nm)';
    'Best-fit sag RMS before PTT (nm)';
    'Best-fit sag PV before PTT (nm)'
    };
Compare_Value = [
    RMS_Direct_PTT_nm;
    PV_Direct_PTT_nm;
    RMS_BestFit_PTT_nm;
    PV_BestFit_PTT_nm;
    RMS;
    PV
    ];
Compare_Table = table(Compare_Item, Compare_Value);
Compare_Table.Compare_Item = char(Compare_Table.Compare_Item);
disp(Compare_Table)

fprintf('\n\n=== 高阶拟合分析（研究用，不建议直接对标 SigFit 前三阶结果）===\n\n');
Zern_Fitting_Results = [RMS; length(zernike_coef); rms(zernike_resid)*1e6; RMS_Zern_PTT_Removed*1e6; PV_Zern_PTT_Removed*1e6; RMS_Zern_PTTPower_Removed*1e6; PV_Zern_PTTPower_Removed*1e6];
XY_Fitting_Results   = [RMS; length(xy_coef);      rms(xy_resid)*1e6;      RMS_XY_PTT_Removed*1e6;   PV_XY_PTT_Removed*1e6;   RMS_XY_PTTPower_Removed*1e6;   PV_XY_PTTPower_Removed*1e6];
Fitting_Iterm = {'Best-fit sag RMS (nm)';'拟合项数';'全拟合残差 RMS (nm)';'去PTT后RMS (nm)';'去PTT后PV (nm)';'去PTTP后RMS (nm)';'去PTTP后PV (nm)'};
Fitting_Matrix = table(Fitting_Iterm, Zern_Fitting_Results, XY_Fitting_Results);
Fitting_Matrix.Fitting_Iterm = char(Fitting_Iterm);
disp(Fitting_Matrix)

%% 10. slope RMS（研究用）
slope_opts_grid.method = 'grid';
slope_opts_grid.grid_nx = 200;
slope_opts_grid.remove_tilt = true;
slope_opts_grid.interp_method = 'natural';

slope_opts_ls.method = 'lsplane';
slope_opts_ls.k_neighbors = 12;
slope_opts_ls.remove_tilt = true;

[sx_rms_g, sy_rms_g, s_mag_rms_g, slopes_grid, info_grid] = compute_slope_rms(Xfit, Yfit, Wz, slope_opts_grid);
[sx_rms_l, sy_rms_l, s_mag_rms_l, slopes_ls, info_ls] = compute_slope_rms(Xfit, Yfit, Wz, slope_opts_ls);

slopeRMS_grid_urad = s_mag_rms_g * 1e6;
slopeRMS_lsplane_urad = s_mag_rms_l * 1e6;

fprintf('\nSlope RMS (grid)    : %.3f urad\n', slopeRMS_grid_urad);
fprintf('Slope RMS (lsplane) : %.3f urad\n', slopeRMS_lsplane_urad);

%% 11. 结果导出
res_t = table(nodeID, Xfit, Yfit, Wz*1e6, zernike_fit*1e6, zernike_resid*1e6, xy_fit*1e6, xy_resid*1e6, ...
    'VariableNames', {'NodeID','X_local_mm','Y_local_mm','BestFitSag_nm','ZernikeFit_nm','ZernikeResid_nm','XYFit_nm','XYResid_nm'});

if isfield(slopes_ls,'sx') && numel(slopes_ls.sx) == numel(nodeID)
    res_t.SlopeX_urad = slopes_ls.sx(:) * 1e6;
    res_t.SlopeY_urad = slopes_ls.sy(:) * 1e6;
    res_t.SlopeMag_urad = sqrt(res_t.SlopeX_urad.^2 + res_t.SlopeY_urad.^2);
end

writetable(res_t, 'mirror_fit_result_sector.csv');
writematrix(zernike_coef, 'zernike_coefs_sector.csv');
writematrix(xy_coef, 'xy_poly_coefs_sector.csv');

summary_table = table( ...
    [1;2;3;4], ...
    {'Direct fit remove PTT'; 'Best-fit plane remove PTT'; 'Zernike remove PTT (analysis)'; 'XY remove PTT (analysis)'}, ...
    [PV_Direct_PTT_nm; PV_BestFit_PTT_nm; PV_Zern_PTT_Removed*1e6; PV_XY_PTT_Removed*1e6], ...
    [RMS_Direct_PTT_nm; RMS_BestFit_PTT_nm; RMS_Zern_PTT_Removed*1e6; RMS_XY_PTT_Removed*1e6], ...
    'VariableNames', {'Index','Title','PV_nm','RMS_nm'});
summary_table.SlopeRMS_grid_urad = repmat(slopeRMS_grid_urad, height(summary_table), 1);
summary_table.SlopeRMS_lsplane_urad = repmat(slopeRMS_lsplane_urad, height(summary_table), 1);
writetable(summary_table, 'plot_pv_rms_summary.csv');

fprintf('\n结果文件已输出：\n');
fprintf('  - mirror_fit_result_sector.csv\n');
fprintf('  - zernike_coefs_sector.csv\n');
fprintf('  - xy_poly_coefs_sector.csv\n');
fprintf('  - plot_pv_rms_summary.csv\n');

%% ========== Zernike 基函数 ==========
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
        fac = norm(Zvec);
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

%% ========== XY 多项式基函数 ==========
function [XYmat, terms] = generateXYPolynomialBasis(x, y, order)
    x = x(:); y = y(:);
    columns = {};
    terms = [];
    for deg = 0:order
        for i = 0:deg
            j = deg - i;
            columns{end+1} = x.^i .* y.^j;
            terms = [terms; i, j];
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

%% ========== slope RMS ==========
function [sx_rms, sy_rms, s_mag_rms, slopes, info] = compute_slope_rms(X, Y, Z_mm, opts)
    % 统一使用 mm；由于 slope 是 dZ/dX，保持一致单位即可。
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'method'), opts.method = 'grid'; end
    if ~isfield(opts,'grid_nx'), opts.grid_nx = 200; end
    if ~isfield(opts,'remove_tilt'), opts.remove_tilt = true; end
    if ~isfield(opts,'k_neighbors'), opts.k_neighbors = 12; end
    if ~isfield(opts,'interp_method'), opts.interp_method = 'natural'; end

    Xc = X(:); Yc = Y(:); Zc = Z_mm(:);
    info = opts;

    if opts.remove_tilt
        A = [Xc, Yc, ones(size(Xc))];
        coeff = A \ Zc;
        Zc = Zc - A * coeff;
        info.removed_plane = coeff;
    else
        info.removed_plane = [0;0;0];
    end

    switch lower(opts.method)
        case 'grid'
            nx = opts.grid_nx;
            xgv = linspace(min(Xc), max(Xc), nx);
            ygv = linspace(min(Yc), max(Yc), nx);
            [Xg, Yg] = meshgrid(xgv, ygv);
            F = scatteredInterpolant(Xc, Yc, Zc, opts.interp_method, 'none');
            Zg = F(Xg, Yg);
            valid_mask = ~isnan(Zg) & ~isinf(Zg);
            if nnz(valid_mask) == 0
                error('插值后无有效网格点，请改用 lsplane 方法或调整参数。');
            end
            dx = xgv(2) - xgv(1);
            dy = ygv(2) - ygv(1);
            [dZdy, dZdx] = gradient(Zg, dy, dx);
            sx_all = dZdx;
            sy_all = dZdy;
            sx = sx_all(valid_mask);
            sy = sy_all(valid_mask);
            s_mag = sqrt(sx.^2 + sy.^2);
            sx_rms = sqrt(mean(sx.^2));
            sy_rms = sqrt(mean(sy.^2));
            s_mag_rms = sqrt(mean(s_mag.^2));
            slopes.sx = sx_all; slopes.sy = sy_all; slopes.mask = valid_mask;
            slopes.Xg = Xg; slopes.Yg = Yg; slopes.Zg = Zg;
            info.method = 'grid';
            info.grid_dx = dx; info.grid_dy = dy;
        case 'lsplane'
            pts = [Xc, Yc];
            M = size(pts,1);
            K = min(opts.k_neighbors, M);
            idxs = knnsearch(pts, pts, 'K', K);
            a_vec = zeros(M,1);
            b_vec = zeros(M,1);
            for i = 1:M
                nb = idxs(i,:);
                Xi = Xc(nb); Yi = Yc(nb); Zi = Zc(nb);
                G = [Xi, Yi, ones(numel(Xi),1)];
                c = G \ Zi;
                a_vec(i) = c(1);
                b_vec(i) = c(2);
            end
            sx = a_vec; sy = b_vec;
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

    info.sx_rms_rad = sx_rms;
    info.sy_rms_rad = sy_rms;
    info.s_mag_rms_rad = s_mag_rms;
    info.sx_rms_urad = sx_rms * 1e6;
    info.sy_rms_urad = sy_rms * 1e6;
    info.s_mag_rms_urad = s_mag_rms * 1e6;
end
