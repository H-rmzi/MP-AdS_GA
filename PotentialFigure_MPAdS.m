%% Reproduce_PotentialFigure_ExactStyle.m
%Fig. 1. paper: The Singly Spinning Myers--Perry--AdS Gravitational Atom:
%  Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths
% MP-AdS, D=7, K=2, l=0 and l=2:  the TWO-PANEL
% analytical effective-potential figure.
% Execute this file as a SCRIPT, from any folder. The PDF, PNG and CSV
% are saved in the current working directory.
clear; clc; close all;

%% PARAMETERS: EXACTLY THOSE USED FOR THE ARTICLE'S STATIC POTENTIAL
D     = 7;
ell   = 50;
mu_s  = 0.1;
h     = 0.810513;          % r_m^(D-3), nonrotating metric parameter
K     = 2;

a     = 0;                 % Figure is STATIC, even though the paper is rotating
styles = [177 85 49; 36 85 138]/255; % orange: l=0, blue: l=2
%styles = [0.5 0 0; 0 0 0.5];
% Dimensional metric and its derivative, in geometric G=c=1 units.
F  = @(r) 1 + r.^2/ell^2 - h./r.^(D-3);
Fp = @(r) 2*r/ell^2 + (D-3)*h./r.^(D-2);

% Dimensional, frequency-independent Schrödinger potential, a=0.
V = @(r,l) F(r).*(mu_s^2 + l*(l+D-3)./r.^2 ...
           + (D-2)*(D-4)*F(r)./(4*r.^2) + (D-2)*Fp(r)./(2*r));

% Exact GLOBAL-AdS scalar frequency for K=2. Do NOT replace by 0.2.
DeltaPlus = (D-1)/2 + sqrt((D-1)^2/4 + (mu_s*ell)^2);
omega_AdS = (DeltaPlus+K)/ell;
E        = omega_AdS^2;

%% HORIZON AND TORTOISE COORDINATE: r_*(r=1.2)=0
rplus = fzero(F,[0.7 1.2]);
fp    = Fp(rplus);
rref  = 1.2;

% Dense, irregular mesh resolves the sharp barrier and distant AdS wall.
r = unique(sort([rplus + logspace(-5,log10(0.2),900), ...
                linspace(rplus+0.2,5,570), ...
                linspace(5,135,2300), ...
                linspace(135,320,240), rref]));

% Horizon-regularized integration:
%   dr_*/dr = 1/F,  1/F = 1/[F'(r_+)(r-r_+)] + regular.
% This avoids direct quadrature of a logarithmic singularity.
g = 1./F(r) - 1./(fp*(r-rplus));
J = cumtrapz(r,g);
Jref = interp1(r,J,rref,'linear');
rstar = (J-Jref) + log((r-rplus)/(rref-rplus))/fp;

%% TURNING POINTS WITH RESPECT TO THE SAME AdS REFERENCE ENERGY
turnR = NaN(2,3);
turnX = NaN(2,3);
for il=1:2
    l = 2*(il-1);
    fvalues = V(r,l) - E;
    crossings = find(fvalues(1:end-1).*fvalues(2:end)<0);
    if numel(crossings) ~= 3
        error('Expected 3 turning points for l=%d; found %d.',l,numel(crossings));
    end
    for k=1:3
        ix = crossings(k);
        turnR(il,k) = fzero(@(rr) V(rr,l)-E,[r(ix) r(ix+1)]);
        turnX(il,k) = interp1(r,rstar,turnR(il,k),'linear');
    end
end

%% FIGURE LAYOUT: same two distinct y-scales as the article
fig=figure('Color','w','Units','centimeters','Position',[2.5 2.5 30 13]);
try
    fig.Renderer='painters';
catch
end
T = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1=nexttile(T,1); hold(ax1,'on');
ax2=nexttile(T,2); hold(ax2,'on');
ax = [ax1 ax2];

handles = gobjects(1,3);
for ip=1:2
    A=ax(ip);
    for il=1:2
        l=2*(il-1);
        pot=V(r,l);
        if ip==1
            mask=(rstar>=-3)&(rstar<=77.5)&(pot>0);
        else
            mask=(rstar>=4.5)&(rstar<=60.5)&(pot>0);
        end
        if il==1
            lineStyle='-';
        else
            lineStyle='--';
        end
        H=plot(A,rstar(mask),pot(mask),'Color',styles(il,:), ...
       'LineStyle',lineStyle,'LineWidth',2.5);
        if ip==1
            handles(il)=H;
        end
        for k=1:3
            xt=turnX(il,k);
            if ip==1 || (xt>=4.5 && xt<=60.5)
                plot(A,xt,E,'o','LineStyle','none', ...
                     'MarkerSize',6.5,'LineWidth',1.55, ...
                     'MarkerEdgeColor',styles(il,:), ...
                     'MarkerFaceColor','w','HandleVisibility','off');
            end
        end
    end
    L=plot(A,[-3 77.5],[E E],'--','Color',[0.19 0.25 0.30],'LineWidth',1.25);
    if ip==1; handles(3)=L; end
    A.FontName='Times New Roman';
    A.FontSize=10;
    A.LineWidth=0.8;
    A.TickDir='out';
    A.Box='off';
    A.XGrid='off'; A.YGrid='on';
    A.GridColor=[0.80 0.84 0.87];
    A.GridAlpha=0.42;
    xlabel(A,'Shifted tortoise coordinate  r_*', ...
        'FontName','Times New Roman','FontSize',11.5);
    if ip==1
        ylabel(A,'V_l(r) = \mu_s^2 v_l(r)', ...
            'FontName','Times New Roman','FontSize',12,'Interpreter','tex');
    else
        ylabel(A,'V_l(r)', ...
            'FontName','Times New Roman','FontSize',12,'Interpreter','tex');
    end
end

%% PANEL (a): horizon barrier AND diverging AdS wall
set(ax1,'YScale','log');
xlim(ax1,[-3,77.5]); ylim(ax1,[2e-4,24]);
title(ax1,'(a) Near-horizon barrier and AdS wall', ...
      'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold');
text(ax1,5.9,0.95,'Barrier','Color',[.45 .54 .60], ...
     'FontName','Times New Roman','FontSize',12.5);
text(ax1,62.6,1.32,'AdS wall','Color',[.45 .54 .60], ...
     'FontName','Times New Roman','FontSize',12.5);

%% PANEL (b): resolved outer trapping intervals for both l values
set(ax2,'YScale','linear');
xlim(ax2,[4.5,60.5]); ylim(ax2,[0.02,0.095]);
yticks(ax2,0.02:0.02:0.08);
title(ax2,'(b) Outer trapping intervals (magnified)', ...
      'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold');
text(ax2,52.3,0.085,'(\omega_{K=2}^{AdS})^2 = 0.04692', ...
    'FontName','Times New Roman','FontSize',10.5, ...
    'Interpreter','tex','Color',[.28 .35 .43]);

% Two horizontal colored spans under the wells show the classically
% oscillatory intervals (between the outer two turning points).
barYs=[0.0285,0.0247];
for il=1:2
    xa=turnX(il,2); xb=turnX(il,3); yb=barYs(il);
    plot(ax2,[xa xb],[yb yb],'-','Color',styles(il,:), ...
         'LineWidth',1.7,'HandleVisibility','off');
    plot(ax2,[xa xa],[yb-0.0010 yb+0.0010],'-', ...
         'Color',styles(il,:),'LineWidth',1.3,'HandleVisibility','off');
    plot(ax2,[xb xb],[yb-0.0010 yb+0.0010],'-', ...
         'Color',styles(il,:),'LineWidth',1.3,'HandleVisibility','off');
end

%% SHARED TOP LEGEND: avoids duplicating the legend in two panels
lg = legend(ax1,handles, ...
    {'l = 0   (N = 1)','l = 2   (N = 0)', ...
     '(\omega_{K=2}^{AdS})^2'}, ...
    'Orientation','horizontal','Box','off', ...
    'FontName','Times New Roman','FontSize',11.5, ...
    'Interpreter','tex');
try
    lg.Layout.Tile='north';
catch
    lg.Location='northoutside';
end

% Make all labels / tick numbers and legend use Times New Roman.
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');

%% OUTPUT FILES (vector PDF and high-resolution PNG)
pdfname='MPAdS_SectionIII_Potential_TwoModes.pdf';
pngname='MPAdS_SectionIII_Potential_TwoModes.png';
try
    exportgraphics(fig,pdfname,'ContentType','vector');
    exportgraphics(fig,pngname,'Resolution',450);
catch
    set(fig,'PaperPositionMode','auto');
    print(fig,pdfname,'-dpdf','-vector');
    print(fig,pngname,'-dpng','-r450');
end

%% TURNING-POINT TABLE / NUMERIC REPRODUCIBILITY CHECK
lcol=repelem([0;2],3);
idcol=repmat((1:3)',2,1);
rcol=reshape(turnR.',[],1);
xcol=reshape(turnX.',[],1);
Ecol=repmat(E,6,1);
turnTable=table(lcol,idcol,rcol,xcol,Ecol, ...
 'VariableNames',{'l','turning_index','r_phys','rstar_shifted','V_equals_omega2'});
writetable(turnTable,'Potential_TurningPoints_MATLAB.csv');

disp('Turning points (comparison target: article figure):');
disp(turnTable);
fprintf('r_+ = %.12g; omega_AdS = %.12g; omega_AdS^2 = %.12g\n',rplus,omega_AdS,E);
fprintf('The paper figure uses r_*(r=1.2)=0 and a=0.\n');
fprintf('Saved %s, %s and Potential_TurningPoints_MATLAB.csv\n',pdfname,pngname);
