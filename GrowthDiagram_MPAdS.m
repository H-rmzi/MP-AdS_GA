%% GrowthDiagram_MPAdS.m
%  Fig. 3. 
%  paper: The Singly Spinning Myers--Perry--AdS Gravitational Atom:
%  Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths
% Compare absolute magnitudes on a log axis (panel a) and signs on a linear
% axis for l=2 (panel b). No QBS radial numerical eigensolver is used.
% All frequencies are dimensional in geometrized units (G=c=1).

clear; clc; close all;
% ----------------------- Physical parameters -----------------------
D = 7; q = D-3;
ell = 50; mu_s = 0.1; a = 0.1; h = 0.810513;
DeltaPlus = (D-1)/2 + sqrt(((D-1)/2)^2 + (mu_s*ell)^2);
w_ads = (DeltaPlus+2)/ell;
fun = @(r) (r.^2+a^2).*(1+r.^2/ell^2) - h./r.^(D-5);
rplus = fzero(fun,[0.5,1.5]);
OmegaH = a*(1+rplus^2/ell^2)/(rplus^2+a^2);

dN10 = -h/(12*ell^5)*DeltaPlus*(DeltaPlus-1)*(5*DeltaPlus+26);
dN02 = -h/(24*ell^5)*DeltaPlus*(DeltaPlus-1)*(DeltaPlus-2);
dLTperM = -dN02*(a/ell);

N = [1 0 0 0 0 0];
l = [0 2 2 2 2 2];
m = [0 -2 -1 0 1 2];
omegaR = zeros(1,6); omegaI = zeros(1,6);

for k = 1:numel(N)
    Nk = N(k); lk = l(k); mk = m(k);
    u = lk + q/2; v = DeltaPlus - (q+2)/2;
    normNl = ell^(2*lk+q+2)/2 * gamma(Nk+1)*gamma(u+1)^2*gamma(Nk+v+1) / ...
             (gamma(Nk+u+1)*(2*Nk+u+v+1)*gamma(Nk+u+v+1));
    lambda = lk/q;
    GammaNl = rplus^(2*lk+q+1) / (2*w_ads*normNl) * ...
              (gamma(1+lambda)^4/gamma(1+2*lambda)^2);
    if lk == 0
        omegaR(k) = w_ads + dN10;
    else
        omegaR(k) = w_ads + dN02 + mk*dLTperM;
    end
    omegaI(k) = -GammaNl*(omegaR(k)-mk*OmegaH);
end
assert(omegaI(6)>0 && all(omegaI(1:5)<0));

% ----------------------- Journal-style colors -----------------------
%Cblue = [36 85 138]/255;     % paper blue (decay)
%Cred  = [177 85 49]/255;     % paper rust red (growth)
Cgrid = [0.84 0.86 0.89];
Czero = [0.46 0.55 0.66];
Cblue = [0 0 0.6];     
Cred  = [0.8 0 0];
% Explicit normalized axes placements keep the legend UNDER panel (a),
% and prevent tiledlayout from moving plots when labels are changed.
fig = figure('Color','w','Units','centimeters','Position',[2 2 28.0 13.1]);
ax1 = axes('Parent',fig,'Units','normalized','Position',[0.105 0.240 0.385 0.650]);
ax2 = axes('Parent',fig,'Units','normalized','Position',[0.615 0.240 0.355 0.650]);

% ----------------------- Panel (a): magnitudes -----------------------
hold(ax1,'on');
ypos = 5:-1:0;  % top to bottom (six physical modes)
labels = {'$(1,0,0)$', '$l=2,\,m=-2$', '$l=2,\,m=-1$', ...
          '$l=2,\,m=0$', '$l=2,\,m=+1$', '$l=2,\,m=+2$'};

for k = 1:6
    xx = abs(omegaI(k)); yy = ypos(k);
    if omegaI(k)>0
        cc = Cred;  face = Cred;
    else
        cc = Cblue; face = 'w';
    end
    plot(ax1, xx, yy, 'o', 'MarkerSize',7.8, 'LineWidth',1.5, ...
         'MarkerEdgeColor',cc, 'MarkerFaceColor',face);
end

% Tick VALUES must be ascending even if modes read from top to bottom.
set(ax1,'XScale','log','XLim',[5e-17 2e-7],'YLim',[-0.65 5.65], ...
    'YTick',0:5,'YTickLabel',labels(end:-1:1), ...
    'TickLabelInterpreter','latex', ...
    'XTick',[1e-16 1e-14 1e-12 1e-10 1e-8], ...
    'FontName','Times New Roman','FontSize',11.0, ...
    'Box','off','LineWidth',1.0,'XMinorTick','off','YMinorTick','off', ...
    'XMinorGrid','off','YMinorGrid','off');
ax1.XGrid = 'on'; ax1.YGrid = 'on';
ax1.GridColor = Cgrid; ax1.GridAlpha = 0.4;
ax1.Layer = 'bottom';

% Thin reference separator between the (1,0,0) and l=2 rows.
plot(ax1,[5e-17 2e-7],[4.5 4.5],'--', ...
    'Color',Czero,'LineWidth',0.75,'HandleVisibility','off');

xlabel(ax1,'$|\mathrm{Im}\,\omega|$  (logarithmic scale)', ...
    'Interpreter','latex','FontSize',12.2);
title(ax1,'(a) Growth/decay magnitude hierarchy', ...
    'FontName','Times New Roman','FontWeight','bold','FontSize',13.2);

% Legend OUTSIDE axes, below the x-axis label (as in the manuscript).
hDecay  = plot(ax1,nan,nan,'o-','Color',Cblue,'MarkerFaceColor','w', ...
    'MarkerSize',7.5,'LineWidth',1.7);
hGrowth = plot(ax1,nan,nan,'o-','Color',Cred,'MarkerFaceColor',Cred, ...
    'MarkerSize',7.5,'LineWidth',1.7);
lg = legend(ax1,[hDecay hGrowth], ...
    {'$\mathrm{Im}\,\omega<0$ (decay)', '$\mathrm{Im}\,\omega>0$ (growth)'}, ...
    'Interpreter','latex','Orientation','horizontal','Box','off', ...
    'FontSize',10.5,'AutoUpdate','off');
lg.Units = 'normalized';
lg.Position = [0.112 0.025 0.380 0.056];

% ----------------------- Panel (b): signed rates -----------------------
hold(ax2,'on');

% Each dotted guide is tied to its OWN analytical level.  In particular,
% negative rates descend from Im(omega)=0, while a positive rate rises
% from zero.  Draw guides first so bars and markers remain in front.
for k = 2:6
    xx = m(k);
    yy = 1e14 * omegaI(k);  % same scaled value as the level marker
    plot(ax2,[xx xx],[0 yy],':', ...
        'Color',[0.65 0.73 0.82],'LineWidth',1.0, ...
        'HandleVisibility','off');
end

% Horizontal line corresponding to Im(omega)=0.
yline(ax2,0,'--','Color',Czero,'LineWidth',1.1);

% Exact signed levels, plotted ON TOP of the guide segments.
for k = 2:6
    xx = m(k); yy = 1e14*omegaI(k);
    if omegaI(k)>0
        cc = Cred; face = Cred;
    else
        cc = Cblue; face = 'w';
    end
    plot(ax2,[xx-0.27 xx+0.27],[yy yy],'-', ...
        'Color',cc,'LineWidth',3.0);
    plot(ax2,xx,yy,'o','MarkerSize',8.0,'LineWidth',1.5, ...
        'MarkerEdgeColor',cc,'MarkerFaceColor',face);
end

set(ax2,'XLim',[-2.5 2.6],'YLim',[-1.86 0.41], ...
    'XTick',-2:2,'YTick',[-1.6 -1.2 -0.8 -0.4 0 0.4], ...
    'FontName','Times New Roman','FontSize',11.0, ...
    'Box','off','LineWidth',1.0);
ax2.Layer = 'bottom';

xlabel(ax2,'Azimuthal quantum number $m$', ...
    'Interpreter','latex','FontSize',12.2);
ylabel(ax2,'$10^{14}\,\mathrm{Im}\,\omega_{02m}$', ...
    'Interpreter','latex','FontSize',13.0);
title(ax2,'(b) Signed $l=2$ rates (magnified)', ...
    'Interpreter','latex','FontWeight','bold','FontSize',13.2);
text(ax2,-2.40,0.095,'$\mathrm{Im}\,\omega=0$', ...
    'Interpreter','latex','Color',[0.39 0.48 0.58],'FontSize',11.0);

% Annotation: both text and arrow use AXIS DATA coordinates, not figure
% coordinates. Thus the arrow always points near the growing m=+2 point.
y_growth = 1e14*omegaI(6);
text(ax2,-0.55,0.285, ...
    sprintf('$m=+2: +%.4g\\times10^{-16}$',omegaI(6)*1e16), ...
    'Interpreter','latex','FontName','Times New Roman', ...
    'Color',Cred,'FontSize',11.6);
arrow_x0=1.33; arrow_y0=0.200;
arrow_x1=1.95; arrow_y1=y_growth+0.038;
quiver(ax2,arrow_x0,arrow_y0,arrow_x1-arrow_x0,arrow_y1-arrow_y0,0, ...
    'Color',Cred,'LineWidth',1.2,'MaxHeadSize',0.50,'AutoScale','off', ...
    'HandleVisibility','off');

% The data marker is kept on top of the pointing arrow.
plot(ax2,2,y_growth,'o','MarkerSize',8.0,'LineWidth',1.5, ...
    'MarkerEdgeColor',Cred,'MarkerFaceColor',Cred);

% Reassert explicit plot locations if the legend altered the layout.
ax1.Units='normalized'; ax1.Position=[0.105 0.240 0.385 0.650];
ax2.Units='normalized'; ax2.Position=[0.615 0.240 0.355 0.650];

drawnow;

exportgraphics(fig,'fig_Growth_Decay.pdf','ContentType','vector');
exportgraphics(fig,'fig_Growth_Decay.png','Resolution',450);
GrowthTable = table(N(:),l(:),m(:),omegaR(:),omegaI(:), ...
    'VariableNames',{'N','l','m','omega_real','omega_imag'});
writetable(GrowthTable,'GrowthDiagram_CorrectAnalytic_Data_MATLAB.csv');

fprintf('rplus = %.12g, OmegaH = %.12g\n',rplus,OmegaH);
disp(GrowthTable);

