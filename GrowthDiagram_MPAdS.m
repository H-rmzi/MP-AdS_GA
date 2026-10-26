%% GrowthDiagram_MPAdS.m
%  Fig. 3. paper: The Singly Spinning Myers--Perry--AdS Gravitational Atom:
%  Fine Structure, Rotational Hyperfine Splitting, and Superradiant Widths
% Compare absolute magnitudes on a log axis (panel a) and signs on a linear
% axis for l=2 (panel b). No QBS radial numerical eigensolver is used.
% All frequencies are dimensional in geometrized units (G=c=1).
clear; clc; close all;
D=7; q=D-3; ell=50; mu_s=0.1; a=0.1; h=0.810513;
DeltaPlus=(D-1)/2 + sqrt(((D-1)/2)^2 + (mu_s*ell)^2);
w_ads=(DeltaPlus+2)/ell;
fun=@(r) (r.^2+a^2).*(1+r.^2/ell^2)-h./r.^(D-5);
rplus=fzero(fun,[0.5,1.5]);
% Horizon angular velocity relative to a nonrotating AdS-boundary frame.
OmegaH=a*(1+rplus^2/ell^2)/(rplus^2+a^2);
dN10=-h/(12*ell^5)*DeltaPlus*(DeltaPlus-1)*(5*DeltaPlus+26);
dN02=-h/(24*ell^5)*DeltaPlus*(DeltaPlus-1)*(DeltaPlus-2);
dLTperM=-dN02*(a/ell);
N=[1 0 0 0 0 0]; l=[0 2 2 2 2 2]; m=[0 -2 -1 0 1 2];
omegaR=zeros(1,6);omegaI=zeros(1,6);
for k=1:numel(N)
    Nk=N(k);lk=l(k);mk=m(k);
    u=lk+q/2;v=DeltaPlus-(q+2)/2;
    normNl = ell^(2*lk+q+2)/2 * gamma(Nk+1)*gamma(u+1)^2*gamma(Nk+v+1) / ...
          (gamma(Nk+u+1)*(2*Nk+u+v+1)*gamma(Nk+u+v+1));
    lambda=lk/q;
    GammaNl = rplus^(2*lk+q+1) / (2*w_ads*normNl) * ...
              (gamma(1+lambda)^4/gamma(1+2*lambda)^2);
    if lk==0
        omegaR(k)=w_ads+dN10;
    else
        omegaR(k)=w_ads+dN02+mk*dLTperM;
    end
    omegaI(k)=-GammaNl*(omegaR(k)-mk*OmegaH);
end
assert(omegaI(6)>0 && all(omegaI(1:5)<0));
Cblue=[36 85 138]/255; Cred=[177 85 49]/255;
fig=figure('Color','w','Units','centimeters','Position',[2 2 29 12.5]);
t=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1=nexttile(t,1); hold(ax1,'on');
ypos=5:-1:0;
labels = {'(1,0,0)','l=2,m=-2','l=2,m=-1','l=2,m=0','l=2,m=+1','l=2,m=+2'};
for k=1:6
    x=abs(omegaI(k)); yy=ypos(k);
    if omegaI(k)>0; cc=Cred;else;cc=Cblue;end
    plot(ax1,[x/10^0.11 x*10^0.11],[yy yy],'-','LineWidth',3,'Color',cc);
    plot(ax1,x,yy,'o','MarkerSize',6.5,'MarkerEdgeColor',cc,...
         'MarkerFaceColor','w','LineWidth',1.3);
end
[yt_sorted, idx] = sort(ypos);
set(ax1,'XScale','log','XLim',[5e-17 2.2e-7], 'YLim',[-.65 5.65]);
%set(ax1,'YTick',ypos,'YTickLabel',{'(1,0,0)','l=2,m=-2','l=2,m=-1','l=2,m=0','l=2,m=+1','l=2,m=+2'});
set(ax1,'YTick',yt_sorted,'YTickLabel',labels(idx));
xticks(ax1,[1e-16 1e-14 1e-12 1e-10 1e-8]);grid(ax1,'on');
ax1.YGrid='off'; ax1.Box='off';
xlabel(ax1,'|Im(\omega)| (logarithmic)');
title(ax1,'(a) Growth/decay magnitude hierarchy');

ax2=nexttile(t,2);hold(ax2,'on');
for k=2:6
    xx=m(k); yy=1e14*omegaI(k);
    if omegaI(k)>0; cc=Cred;else;cc=Cblue;end
    plot(ax2,[xx-.28 xx+.28],[yy yy],'-','LineWidth',3,'Color',cc);
    plot(ax2,xx,yy,'o','MarkerSize',6.5,'MarkerEdgeColor',cc,...
         'MarkerFaceColor','w','LineWidth',1.3);
end
yline(ax2,0,'--','Color',[.48 .54 .60]);
set(ax2,'XLim',[-2.5 2.6],'YLim',[-1.89 .41],'XTick',-2:2);
ax2.Box='off'; xlabel(ax2,'Azimuthal quantum number m');
ylabel(ax2,'10^{14} Im(\omega_{02m})');
title(ax2,'(b) Signed l=2 rates (magnified)');
text(ax2,1.15,.27,sprintf('m=+2: +%.4g',omegaI(6)),'Color',Cred,'FontSize',10);

% Times New Roman on systems (e.g. Windows) where it is installed.
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
set(findall(gcf,'Type','axes'),'TickLabelInterpreter','tex');
exportgraphics(fig,'fig_Growth_Decay.pdf','ContentType','vector');
exportgraphics(fig,'fig_Growth_Decay.png','Resolution',450);
GrowthTable=table(N(:),l(:),m(:),omegaR(:),omegaI(:), ...
    'VariableNames',{'N','l','m','omega_real','omega_imag'});
writetable(GrowthTable,'GrowthDiagram_CorrectAnalytic_Data_MATLAB.csv');
fprintf('rplus = %.12g, OmegaH = %.12g\n',rplus,OmegaH);
disp(GrowthTable);
