%% LevelDiagram_MPAdS.m
% Fig. 2. Paper: The Singly Spinning Myers--Perry--AdS Gravitational Atom: 
% Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths
% Uses analytic D=7 (K=2) real-frequency corrections.
% Horizontal segments represent spectral levels. Two vertical scales are
% retained; these are genuine frequency differences.


clear; close all; clc;
D = 7; ell = 50; mu_s = 0.1; a = 0.1;
h = 0.810513;  % h = r_m^4, fixed as in the manuscript's D=7 example
K = 2;
assert(D == 7 && K == 2, 'Closed-form formulas here require D=7, K=2.');

Delta_plus = (D-1)/2 + sqrt(((D-1)/2)^2 + (mu_s*ell)^2);
omega_AdS = (Delta_plus + K)/ell;

% Exact-AdS-wavefunction first-order Newtonian corrections, dimensional omega
dw_10_N = -h/(12*ell^5)*Delta_plus*(Delta_plus-1)*(5*Delta_plus+26);
dw_02_N = -h/(24*ell^5)*Delta_plus*(Delta_plus-1)*(Delta_plus-2);
Delta_FS = dw_02_N - dw_10_N;

% Physical LT correction in the frame nonrotating at the AdS boundary.
% Eq. (D7 dimensionless spectrum): dw_02m = dw_02_N*(1-m*a/ell)
LT_per_m = -dw_02_N*a/ell;
m = (-2:2)';
g = [1;4;10;4;1];
omega_R_02m = omega_AdS + dw_02_N + m*LT_per_m;
omega_R_10  = omega_AdS + dw_10_N;

fprintf('omega_K^AdS = %.15g\n',omega_AdS);
fprintf('dw_10^N = %.12e, dw_02^N = %.12e\n',dw_10_N,dw_02_N);
fprintf('Delta_FS = %.12e, LT per m = %.12e\n',Delta_FS,LT_per_m);

T = table([1;zeros(5,1)], [0;2*ones(5,1)], [0;m], [1;g], ...
    [omega_R_10;omega_R_02m], 'VariableNames', ...
    {'N','l','m','multiplicity','omega_real'});
writetable(T,'MPAdS_Figure2_CorrectAnalyticalData.csv');

figure('Color','w','Units','centimeters','Position',[2 2 30 12]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
blue=[0.14 0.33 0.54]; copper=[0.69 0.33 0.19]; ink=[0.16 0.22 0.28];

% Panel (a) - schematic horizontal x-positions, physically scaled y-values
ax1=nexttile; hold(ax1,'on');
y0=0; y10=1e6*dw_10_N; y02=1e6*dw_02_N;
yline(ax1,0,':','Color',[0.76 0.79 0.82]);
plot(ax1,[.13 .80],[0 0],'-','Color',ink,'LineWidth',3.4);
plot(ax1,[.80 1.44],[0 y10],'--','Color',[0.60 0.66 0.72],'LineWidth',.9);
plot(ax1,[.80 2.73],[0 y02],'--','Color',[0.60 0.66 0.72],'LineWidth',.9);
plot(ax1,[1.44 2.02],[y10 y10],'-','Color',copper,'LineWidth',3.4);
plot(ax1,[2.73 3.30],[y02 y02],'-','Color',blue,'LineWidth',3.4);
plot(ax1,1.73,y10,'o','Color',copper,'MarkerFaceColor','w','MarkerSize',5);
plot(ax1,3.015,y02,'o','Color',blue,'MarkerFaceColor','w','MarkerSize',5);
text(ax1,.47,.08,'$\mathcal{K}=2$','Interpreter','latex','HorizontalAlignment','center');
text(ax1,1.72,y10-.10,'$(N,l)=(1,0)$','Interpreter','latex', ...
    'Color',copper,'HorizontalAlignment','center');
text(ax1,3.015,y02-.13,'$(N,l)=(0,2)$','Interpreter','latex', ...
    'Color',blue,'HorizontalAlignment','center');
plot(ax1,[2.43 2.43],[y10 y02],'-','Color',[.40 .46 .52],'LineWidth',1);
text(ax1,2.27,(y10+y02)/2,'$\Delta\omega_{\rm FS}$', ...
    'Rotation',90,'Interpreter','latex','HorizontalAlignment','center');
text(ax1,2.42,-1.21,'$9.975\times10^{-7}$','Interpreter','latex', ...
    'HorizontalAlignment','center');
ax1.XLim=[-.04 3.52]; ax1.YLim=[-1.30 .235];
ax1.XTick=[.46 1.72 3.02];
ax1.XTickLabel={'Pure AdS','$l=0$','$l=2$'};
ylabel(ax1,'$10^{6}\,[\mathrm{Re}\,\omega-\omega_{\mathcal{K}}^{\mathrm{AdS}}]$', ...
    'Interpreter','latex');
title(ax1,'(a) Pure AdS and Newtonian fine structure');
ax1.TickLabelInterpreter='latex'; ax1.Box='off';

% Panel (b) - real hyperfine enlargement, relative to l=2 centroid
ax2=nexttile; hold(ax2,'on');
yline(ax2,0,':','Color',[0.76 0.79 0.82]);
for k=1:numel(m)
    y=1e10*m(k)*LT_per_m;
    plot(ax2,[.07 .45],[y y],'-','Color',blue,'LineWidth',3);
    plot(ax2,.27,y,'o','Color',blue,'MarkerFaceColor','w','MarkerSize',5);
    text(ax2,.50,y,sprintf('$m=%+d\; (g=%d)$',m(k),g(k)), ...
        'Interpreter','latex','VerticalAlignment','middle','FontSize',11);
end
ax2.XLim=[0 1.22]; ax2.YLim=[-2.68 2.72];ax2.XTick=[];
ax2.YTick=[-2.042 -1.021 0 1.021 2.042];
ylabel(ax2,'$10^{10}\,\delta\omega_{02m}^{\mathrm{LT}}$', ...
    'Interpreter','latex');
title(ax2,'(b) Hyperfine sublevels (magnified)');
text(ax2,.5,-2.50, ...
    '$\omega_{02}^{\mathrm{center}}=\omega_{\mathcal{K}}^{\mathrm{AdS}}+\delta\omega_{02}^{\mathrm{N}}$', ...
    'Interpreter','latex','HorizontalAlignment','center','FontSize',10);
ax2.TickLabelInterpreter='latex'; ax2.Box='off';


% Times New Roman on systems (e.g. Windows) where it is installed.
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
set(findall(gcf,'Type','axes'),'TickLabelInterpreter','tex');
exportgraphics(gcf,'MPAdS_Figure2_CorrectAnalytic.pdf','ContentType','vector');
exportgraphics(gcf,'MPAdS_Figure2_CorrectAnalytic.png','Resolution',500);
