function output = MPAdS_Analytical_Supplement(runNodeFree, beta_muL)
% MPADS_ANALYTICAL_SUPPLEMENT One-file for MP--AdS scalar modes.
%
%   MPAdS_Analytical_Supplement
%       Computes static Q_out, (A_N,B_N) for N=0,1 and the independent
%       node-free radial check. Writes output tables to a results folder.
%   MPAdS_Analytical_Supplement(true, 5)
%       Runs all reference tests with beta=mu_s*ell=5.
%   R = MPAdS_Analytical_Supplement(...)
%       Returns the computed coefficients in a struct.
% Conventions: D=7, l=2; H=h/ell^4, xi=m*a/ell, Omega=ell*omega_stat.
% S_N(H) = d(Re Omega)/d xi|_{xi=0}
%  = H*Lambda_N + H^2*[A_N*log(r_plus_static/ell)+B_N]+...
% Written as a *single* .m file with local/nested functions.
% This analytical coefficients; it is NOT a full QBS solver.

if nargin < 1 || isempty(runNodeFree), runNodeFree=true; end
if nargin < 2 || isempty(beta_muL), beta_muL=5; end
validateattributes(runNodeFree,{'logical','numeric'},{'scalar'});
runNodeFree=logical(runNodeFree);
validateattributes(beta_muL,{'numeric'},{'scalar','real','finite','nonnegative'});

root=fileparts(mfilename('fullpath'));
outDir=fullfile(root,'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('MP--AdS analytical coefficient supplement (one MATLAB file)\n');
fprintf('D=7, l=2; beta=mu_s*ell=%.12g; N=0,1\n', beta_muL);

q=mpads_static_qout(beta_muL);
fprintf('Q_out (closed):    %.12f\n',q.Q_closed);
fprintf('Q_out (beta form): %.12f\n',q.Q_from_beta);
assert(q.absolute_difference<5e-9,'Static Q_out formulas disagree.');

c0=mpads_mixed_coefficients(beta_muL,0);
c1=mpads_mixed_coefficients(beta_muL,1);
T=table([c0.N;c1.N],[c0.Omega0;c1.Omega0], ...
    [c0.Omega1;c1.Omega1],[c0.Lambda;c1.Lambda], ...
    [c0.norm;c1.norm],[c0.I_FP;c1.I_FP], ...
    [c0.A;c1.A],[c0.B;c1.B], ...
    'VariableNames',{'N','Omega0','Omega1','Lambda', ...
    'AdSNorm','FarFinitePart','A_log','B_finite'});
disp(T);
writetable(T,fullfile(outDir,'analytical_coefficients.csv'));

if abs(beta_muL-5)<1e-12
    referenceA=[-2095.821341072154;-19683.481322551448];
    referenceB=[815.615518133988;5373.80684372628];
    referenceI=[0.304182657710;1.850688889];
    assert(all(abs(T.A_log-referenceA)<1e-3),'A_N reference check failed.');
    assert(all(abs(T.B_finite-referenceB)<0.05),'B_N reference check failed.');
    assert(all(abs(T.FarFinitePart-referenceI)<1e-4), ...
        'Far finite-part reference check failed.');
    fprintf('PASS: static and mixed-coefficient reference tests.\n');
else
    fprintf('INFO: reference assertions are defined for beta=5 only.\n');
end

nodeTable=table();
if runNodeFree
    fprintf('Running independent node-free radial check...\n');
    nodeTable=mpads_nodefree_finitepart(beta_muL,[.05;.03;.02;.015;.01]);
    disp(nodeTable);
    writetable(nodeTable,fullfile(outDir,'nodefree_convergence_matlab.csv'));
    if abs(beta_muL-5)<1e-12
        referenceNode=1.850688889;
        last=nodeTable.EpsilonSquaredCorrected(end);
        fprintf('At eps=%.4f: I_FP=%.10f, reference=%.10f\n', ...
            nodeTable.epsilon(end),last,referenceNode);
        assert(abs(last-referenceNode)<0.002, ...
            'Node-free check did not converge to the specified tolerance.');
        fprintf('PASS: independent node-free finite-part check.\n');
    end
else
    fprintf('Node-free check skipped (runNodeFree=false).\n');
end

fid=fopen(fullfile(outDir,'run_summary.txt'),'w');
if fid<0, error('Unable to write run_summary.txt.'); end
fprintf(fid,'MATLAB release: %s\n', version);
fprintf(fid,'Source: MPAdS_Analytical_Supplement.m\n');
fprintf(fid,'beta=%.15g, D=7, l=2\n',beta_muL);
fprintf(fid,'Q_out closed=%.15g, reconstructed=%.15g\n', ...
    q.Q_closed,q.Q_from_beta);
for k=1:2
    fprintf(fid,'N=%d A=%.15g B=%.15g I_FP=%.15g\n', ...
        T.N(k),T.A_log(k),T.B_finite(k),T.FarFinitePart(k));
end
fprintf(fid,'Node-free check executed: %d\n',runNodeFree);
if runNodeFree
    fprintf(fid,'Node-free I_FP (smallest epsilon): %.15g\n', ...
        nodeTable.EpsilonSquaredCorrected(end));
end
fclose(fid);
fprintf('Saved results in %s\n',outDir);
output=struct('static',q,'coefficients',T,'nodefree',nodeTable, ...
    'beta_muL',beta_muL,'nodefree_executed',runNodeFree);
end


%% Local helper: mpads_static_qout.m
function out = mpads_static_qout(beta_muL)
%MPADS_STATIC_QOUT Static D=7, (N,l)=(0,2) finite-part coefficient.
% beta_muL = mu_s * ell. Requires MATLAB base functions only.
% Both the closed digamma formula and independent beta-integral
validateattributes(beta_muL, {'numeric'}, {'scalar','real','finite','nonnegative'});
d = 3 + sqrt(9 + beta_muL^2);
Eg = -psi(1); % Euler--Mascheroni constant
out.beta_muL = beta_muL;
out.Delta = d;
out.Q_closed = 19*d^2/576 + 71*d/144 + 31*d/(96*(d+1)) ...
             - (psi(d+1)+Eg)/2;

% An independent beta-function reconstruction of the finite part.
a = (d^2+8*d+6)/6;
b = d*(d+2)*(d+3)/24;
c = d*(d-1)*(d+2)/24;
% Q(y)=(1-y)^3*(-1+a*x^2+b*x^4+c*x^6), y=x^2/(1+x^2).
% Ascending polynomial order for Q(y).
q = [-1, 3+a, -3-2*a+b, 1+a-b+c];
qq = conv(q,q);
finite_integral = -(psi(d-1)+Eg)/2;
for k=1:6
    finite_integral = finite_integral + 0.5*qq(k+1)*exp(betaln(k,d-1));
end
Omega0 = d+2;
Omega1 = -d*(d-1)*(d-2)/24;
w2I = 0.5*exp(betaln(1,d+4));
w1I = 0.5*exp(betaln(3,d+1));
w0I = 0.5*exp(betaln(5,d-2));
out.Q_from_beta = 0.75 + finite_integral + Omega0^2*w2I ...
    + 2*Omega0*Omega1*w1I + Omega1^2*w0I;
out.absolute_difference = abs(out.Q_closed - out.Q_from_beta);
end


%% Local helper: mpads_mixed_coefficients.m
function out = mpads_mixed_coefficients(beta_muL, N)
%MPADS_MIXED_COEFFICIENTS D=7, l=2 spin-odd mixed coefficients N=0,1.
% Conventions: H=h/ell^4, xi=m*a/ell, Omega=ell*omega_stat,
% d Re(Omega)/d xi|_0 = H*Lambda + H^2*(A*log(rplus/ell)+B)+...
% A and B are coefficients in a formal small-hole/low-frequency,
% spin-linear matching expansion (not exact finite-spin QBS frequencies).
validateattributes(beta_muL, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(N, {'numeric'}, {'scalar','integer','>=',0,'<=',1});
d = 3 + sqrt(9+beta_muL^2);
Om = d+2+2*N;
norm0 = 0.5*exp(betaln(5,d-2));
if N==0
    Om1 = -d*(d-1)*(d-2)/24;
    Lambda = d*(d-1)*(d-2)/24;
    norm = norm0;
else
    Om1 = -d*(d-1)*(7*d+46)/120;
    Lambda = d*(d-1)*(11*d+38)/120;
    norm = norm0*(d-2)/(5*(d+4));
end
I_FP = mpads_mixed_far_finitepart(d,N,Om,Om1,Lambda);
zeta3 = 1.20205690315959428539973816151144999;
Knear = -(log(2)+3/8-7*zeta3/16)/8;
A = -1/(2*norm);
B = -(Om*(8*Knear+1/8)+I_FP)/(2*Om*norm);
out.beta_muL = beta_muL;
out.N = N;
out.Delta = d;
out.Omega0 = Om;
out.Omega1 = Om1;
out.Lambda = Lambda;
out.norm = norm;
out.I_FP = I_FP;
out.Knear = Knear;
out.A = A;
out.B = B;
end


%% Local helper: mpads_mixed_far_finitepart.m
function I_FP = mpads_mixed_far_finitepart(d,N,Om,Om1,Lambda)
%MPADS_MIXED_FAR_FINITEPART MATLAB quadrature in y=x^2/(1+x^2).
% The flux primitives J(y), K(y) follow from polynomial first-order
% inhomogeneous radial equations. Both origin log and internal N=1
% double pole are removed analytically before quadrature.
%
% Definition: regular(y) = original(y)-alpha/(y-y0)^2
%                         + Om*(1-y)^(d-2)/(2*y).
% The double-pole finite part on (0,1) is
% -alpha*(1/y0+1/(1-y0)).
% N=0 uses y0 as an ordinary interval splitting point only.
if N==0
    c=0; y0=0.4; alpha=0;
    j = [-24, 4*d^2+32*d+96, d^3-3*d^2-58*d-120,24*d+48]/24;
    k = -(d+2)*[12,4*d-24,d^2-5*d+12]/24;
else
    c=(d+3)/5;
    y0=1/c;
    j = [-3000, 500*d^2+8000*d+24000, ...
        -25*d^3-2465*d^2-25690*d-55560, ...
        -28*d^4-70*d^3+4486*d^2+32812*d+55440, ...
        7*d^5+81*d^4-49*d^3-4125*d^2-18714*d-25200, ...
        120*d^3+1200*d^2+3960*d+4320]/3000;
    k = -(d+4)*[1500,100*d-4200,55*d^2+55*d+4170, ...
        -44*d^3-174*d^2-106*d-1740, ...
        11*d^4+49*d^3+27*d^2+27*d+270]/3000;
    alpha = y0^(-1)*(1-y0)^(d-2) ...
            *polyval(fliplr(j),y0)*polyval(fliplr(k),y0)/(c*c);
end

opts={'RelTol',2e-10,'AbsTol',2e-11};
if N==0
    integ=quadgk(@regfun,0,y0,opts{:})+quadgk(@regfun,y0,1,opts{:});
    boundary=0;
else
    % Avoid loss of significance in double-pole cancellation by
    % evaluating a small central interval using Gauss--Legendre nodes.
    gap=0.01;
    left=quadgk(@regfun,0,y0-gap,opts{:});
    right=quadgk(@regfun,y0+gap,1,opts{:});
    [z,w]=mpads_gauss_legendre(64);
    middle=gap*sum(w.*regfun(y0+gap*z));
    integ=left+middle+right;
    boundary=-alpha*(1/y0+1/(1-y0));
end
Eg=-psi(1);
I_FP=Om/4+integ+boundary+Om/2*(psi(d-1)+Eg);

% Nested numerical evaluator, vectorized for quadgk.
    function f = regfun(y)
        t=1-y;
        py=1-c*y;
        J=polyval(fliplr(j),y);
        K=polyval(fliplr(k),y);
        jk=t.^(d-2).*J.*K./(y.*py.^2);
        direct=py.^2.*(Om1*Lambda*y.^4.*t.^(d-3) ...
            +Om*Lambda*y.^2.*t.^d ...
            -Om1*y.*t.^d -Om*t.^(d+3)./y);
        f=jk+direct +Om*t.^(d-2)./(2*y);
        if N==1
            f=f-alpha./(y-y0).^2;
        end
    end

end


%% Local helper: mpads_nodefree_finitepart.m
function result = mpads_nodefree_finitepart(beta_muL, epsValues)
%MPADS_NODEFREE_FINITEPART Optional independent x-space N=1 audit.
% Solve the mass-driven first-order inhomogeneous radial ODE for R_M,
% avoiding divisions by the nodal AdS eigenfunction. The result is the
% same mixed mass/spin far-region finite part obtained by Hadamard FP.
% Base MATLAB integration functions only (ode113, integral).
% This is a check of far-region perturbation theory, NOT QBS eigenvalues.
if nargin<1, beta_muL=5; end
if nargin<2, epsValues=[0.05;0.03;0.02;0.015;0.01]; end
epsValues=epsValues(:);
d=3+sqrt(9+beta_muL^2);
Om=d+4;
Om1=-d*(d-1)*(7*d+46)/120;
LT=d*(d-1)*(11*d+38)/120;
T=(d+2)/2+(d+3)/5;
C=Om*(14+8*T-Om^2)/6;
X=12;
raw=zeros(size(epsValues)); corrected=raw; wrErr=raw;
opts=odeset('RelTol',2e-10,'AbsTol',1e-13,'MaxStep',0.035);
for ii=1:numel(epsValues)
    eps0=epsValues(ii);
    jM0=integral(@(xx) origin_flux(xx,1),0,eps0,...
        'RelTol',5e-12,'AbsTol',2e-15);
    jS0=integral(@(xx) origin_flux(xx,2),0,eps0,...
        'RelTol',5e-12,'AbsTol',2e-15);
    [r0,~,~]=radial_mode(eps0);
    p0=eps0^5*(1+eps0^2);
    init=[0; -jM0/(p0*r0); jM0; jS0;0;0];
    [~,Y]=ode113(@radial_rhs,[eps0 X],init,opts);
    yF=Y(end,:).';
    [rX,rpX,~]=radial_mode(X);
    pX=X^5*(1+X^2);
    raw(ii)=Om/4+yF(5)+yF(6)-2*yF(4)*yF(1)/rX-Om*log(eps0);
    corrected(ii)=raw(ii)+C*eps0^2/2;
    wrErr(ii)=pX*(rX*yF(2)-rpX*yF(1))+yF(3);
end
result=table(epsValues,raw,corrected,wrErr,'VariableNames',...
    {'epsilon','RawFinitePart','EpsilonSquaredCorrected','WronskianError'});

    function [R,Rp,Rpp] = radial_mode(x)
        % Exact D=7, (N,l)=(1,2) normalizable pure-AdS eigenfunction.
        g=(d+2)/2;
        c=(d+3)/5;
        f=1+x.^2;
        S=f.^(-g);
        Sp=-2*g*x./f.*S;
        Spp=S.*(-2*g./f+4*g*(g+1)*x.^2./f.^2);
        A=x.^2.*S;
        Ap=2*x.*S+x.^2.*Sp;
        App=2*S+4*x.*Sp+x.^2.*Spp;
        y=x.^2./f;
        yp=2*x./f.^2;
        ypp=2*(1-3*x.^2)./f.^3;
        P=1-c*y;
        Pp=-c*yp;
        Ppp=-c*ypp;
        R=A.*P;
        Rp=Ap.*P+A.*Pp;
        Rpp=App.*P+2*Ap.*Pp+A.*Ppp;
    end

    function [FM,FS,direct,V] = sources(x)
        [R,Rp,Rpp]=radial_mode(x);
        f=1+x.^2;
        w0=x.^5./f;
        w1=x./f.^2;
        v1=1./(x.*f);
        v2=1./(x.^5.*f.^2);
        FM=Rp+x.*Rpp-(Om^2*w1+2*Om*Om1*w0).*R;
        FS=-2*Om*(LT*w0-v1).*R;
        direct=(2*Om1*LT*w0+2*Om*LT*w1-2*Om1*v1-2*Om*v2).*R.^2;
        V=Om^2*w0-beta_muL^2*x.^5-12*x.^3;
    end

    function flux = origin_flux(xx,which)
        [FM,FS]=sources(xx);
        R=radial_mode(xx);
        if which==1
            flux=-R.*FM;
        else
            flux=-R.*FS;
        end
    end

    function rhs = radial_rhs(x,y)
        [FM,FS,di,V]=sources(x);
        R=radial_mode(x);
        p=x^5*(1+x^2);
        pp=5*x^4+7*x^6;
        rhs=[y(2); (FM-pp*y(2)-V*y(1))/p; ...
             -R*FM; -R*FS; -2*FS*y(1); di];
    end
end


%% Local helper: mpads_gauss_legendre.m
function [x,w] = mpads_gauss_legendre(n)
%MPADS_GAUSS_LEGENDRE n-point Gauss--Legendre rule on [-1,1].
validateattributes(n,{'numeric'},{'scalar','integer','>=',2});
j=1:n-1;
beta=j./sqrt(4*j.^2-1);
T=diag(beta,1)+diag(beta,-1);
[V,D]=eig(T);
[x,idx]=sort(diag(D));
w=2*(V(1,idx).').^2;
x=x(:);w=w(:);
end
