clc
clear all
close all

%% Number of Loops

n = 10000;

%% Bounds

LB = [ 0,  2.57,  0,   0,   0];      % e(1), a(1), t(1), p(1), p(2)
UB = [50, 13.6 , 42, inf, inf];

%% Utility Coefficients

UCoef = [- 0.4591,  % Price      (10,000 2006$)
         -36.7700,  % 1/econ     (econ)
          11.2620,  % 1/acc      (s)
           2.4541,  % footprint  (10,000 in^2)
         - 8.0903]; % const
     
%% Cost Coefficients

cCoef = [ 0.7800,  % const
          1.9716,  % exp(-acc)  (s)
          0.0016,  % tech
          0.2250,  % wt         (1,000 lbf)
         -0.0123]; % wt*acc     (1,000 lbf, s)

%% Constraint Coefficients

ACoef = [  0.3783,  % 
          10.7920,  % const
          69.5244,  % exp(-acc)   (s, #)
          12.9897,  % wt          (1,000 lbf)
         - 0.5593,  % wt*acc      (1,000 lbf, s)
         - 0.2605,  % tech
           0.0013]; % acc^2*tech  (s, #)
       
%% Fixed Characteristics

eqcost = 1.2668716427791604; % .859087231474419; % Eq. Cost  (10,000 2006$)

e = [  0      , 23.53   ];
a = [  0      ,  6.33   ];
t = [  0      ,  0      ];
p = [  0      ,  0      ];
w = [  3.285  ,  3.285  ]; % Weight    (1,000 lbf)
f = [ 12.7413 , 12.7413 ]; % Footprint (10,000 in^2)
c = [  0      , eqcost  ];


%% Initial Conditions

e0 = LB(1) + (UB(1)-LB(1)).*rand(n,1);  % Econ (mpg)
a0 = LB(2) + (UB(2)-LB(2)).*rand(n,1);  % Acc  (0-60)
t0 = LB(3) + (UB(3)-LB(3)).*rand(n,1);  % Tech

% Prices start below $100000, as this was the cutoff in the Whitefoot paper
p0 = 10 * rand(n,2);        % Price 2

%% Linear Inequality Constraints

A = [];
B = [];

%% Linear Equality Constraints

Aeq = [];
Beq = [];

%% Options

opt = optimset('fmincon');

opt.Algorithm = 'interior-point';
opt.AlwaysHonorConstraints = 'bounds';
opt.InitBarrierParam = 0.1;
opt.InitTrustRegionRadius = 1; % I believe
opt.TolCon = 10^(-6);
opt.FinDiffType = 'forward';
opt.MaxIter = 10000;
opt.TolFun = 10^(-6);
opt.TolX = 10^(-14);

opt.Hessian = 'bfgs';
opt.MaxProjCGIter = 10;

% User-supplied gradients
opt.DerivativeCheck = 'off';
opt.GradConstr = 'on';
opt.GradObj = 'on';

%% Optimization
i = 1;

global Safety
Safety = 1;

for i = 1:n
    
    i
    
    x0 = [ e0(i,1) , a0(i,1) , t0(i,1) , p0(i,1) , p0(i,2) ];
    
    tic
    
    [X,FVAL,EXITFLAG(i,1),OUTPUT(i,1),LAMBDA(i,1)] ...
        = fmincon( @(x)( OBJ(x,e,a,t,w,f,cCoef,UCoef) ),...
          x0,...
          A,...
          B,...
          Aeq,...
          Beq,...
          LB,...
          UB,...
          @(x)( NONLCON(x,e,a,t,w,f,c,ACoef,UCoef,cCoef) ),...
          opt );
    
    time(i,1) = toc;  
      
    es(i,1) = X(1);     es(i,2) = e(2);
    as(i,1) = X(2);     as(i,2) = a(2);
    ts(i,1) = X(3);     ts(i,2) = t(2);
    ps(i,1:2) = X(4:5);
    
    Us(i,1) = UCoef(1)* ps(i,1)... % 10,000 2006$
            + UCoef(2)/ es(i,1)... % mpg
            + UCoef(3)/ as(i,1)... % s
            + UCoef(4)* f(1)...    % 10,000 in^2
            + UCoef(5);
    
    Us(i,2) = UCoef(1)* ps(i,2)... % 10,000 2006$
            + UCoef(2)/ es(i,2)... % mpg
            + UCoef(3)/ as(i,2)... % s
            + UCoef(4)* f(2)...    % 10,000 in^2
            + UCoef(5);
    
    Ps(i,1) = exp( Us(i,1) ) / ( 1.0 + exp( Us(i,1) ) + exp( Us(i,2) ) );
    Ps(i,2) = exp( Us(i,2) ) / ( 1.0 + exp( Us(i,1) ) + exp( Us(i,2) ) );
    
    cs(i,1) = cCoef(1)...
            + cCoef(2) * exp( -as(i,1) )...
            + cCoef(3) * ts(i,1) ...
            + cCoef(4) * w(1) ...
            + cCoef(5) * w(1) * as(i,1);
    cs(i,2) = c(2);
    
    pr(i,1) = Ps(i,1) * ( ps(i,1) - cs(i,1) );
    
end

%% Optimization (Zeta)
i = 1;

for i = 1:n
    
    i
    
    x0 = [ e0(i,1) , a0(i,1) , t0(i,1) , p0(i,1) , p0(i,2) ];
    
    tic
    
    [XZ,FVALZ,EXITFLAGZ(i,1),OUTPUTZ(i,1),LAMBDAZ(i,1)] ...
        = fmincon( @(x)( OBJ(x,e,a,t,w,f,cCoef,UCoef) ),...
          x0,...
          A,...
          B,...
          Aeq,...
          Beq,...
          LB,...
          UB,...
          @(x)( NONLCONZ(x,e,a,t,w,f,c,ACoef,UCoef,cCoef) ),...
          opt );
    
    timeZ(i,1) = toc;  
      
    esZ(i,1) = XZ(1);     esZ(i,2) = e(2);
    asZ(i,1) = XZ(2);     asZ(i,2) = a(2);
    tsZ(i,1) = XZ(3);     tsZ(i,2) = t(2);
    psZ(i,1:2) = XZ(4:5);
    
    UsZ(i,1) = UCoef(1)* psZ(i,1)... % 10,000 2006$
            + UCoef(2)/ esZ(i,1)...  % mpg
            + UCoef(3)/ asZ(i,1)...  % s
            + UCoef(4)* f(1)...      % 10,000 in^2
            + UCoef(5);
    
    UsZ(i,2) = UCoef(1)* psZ(i,2)... % 10,000 2006$
            + UCoef(2)/ esZ(i,2)...  % mpg
            + UCoef(3)/ asZ(i,2)...  % s
            + UCoef(4)* f(2)...      % 10,000 in^2
            + UCoef(5);
    
    PsZ(i,1) = exp( UsZ(i,1) ) / ( 1.0 + exp( UsZ(i,1) ) + exp( UsZ(i,2) ) );
    PsZ(i,2) = exp( UsZ(i,2) ) / ( 1.0 + exp( UsZ(i,1) ) + exp( UsZ(i,2) ) );
    
    csZ(i,1) = cCoef(1)...
            + cCoef(2) * exp( -asZ(i,1) )...
            + cCoef(3) * tsZ(i,1) ...
            + cCoef(4) * w(1) ...
            + cCoef(5) * w(1) * asZ(i,1);
    csZ(i,2) = c(2);
    
    prZ(i,1) = PsZ(i,1) * ( psZ(i,1) - csZ(i,1) );
    
end

%% Bar Plot

outcome = zeros(1,3); %Success, False success, Failed
Success = zeros(n,1);

for i = 1:n

    if ( EXITFLAG(i) == 1 || EXITFLAG(i) == 2 )
        
        if( Ps(i,1) > 1e-6 && Ps(i,2) > 1e-6 )
        
            outcome(1,1) = outcome(1,1) + 1;
            Success(i) = 1;
                        
        else
            
            outcome(1,2) = outcome(1,2) + 1;
            
        end
        
    else 
        
        outcome(1,3) = outcome(1,3) + 1;
        
    end

end

outcome = outcome/n;

figure(1)

bar(outcome)

title('Outcome Plot')
ylabel('Portion of Runs')
set(gca,'XTickLabel',{'Success', 'False Success', 'Failure'})
set(gca,'YLim',[0,1]);

Sind = find(Success);

%saveas(gcf,'BarPlot.jpg')

%% Bar Plot (Zeta)

outcomeZ = zeros(1,3); %Success, False success, Failed
SuccessZ = zeros(n,1);

for i = 1:n

    if ( EXITFLAGZ(i) == 1 || EXITFLAGZ(i) == 2 )
        
        if( PsZ(i,1) > 1e-6 && PsZ(i,2) > 1e-6 )
        
            outcomeZ(1,1) = outcomeZ(1,1) + 1;
            SuccessZ(i) = 1;
            
        else
            
            outcomeZ(1,2) = outcomeZ(1,2) + 1;
            
        end
        
    else 
        
        outcomeZ(1,3) = outcomeZ(1,3) + 1;
        
    end

end

outcomeZ = outcomeZ/n;

figure(2)

bar(outcomeZ)

title('Outcome Plot (Zeta)')
ylabel('Portion of Runs')
set(gca,'XTickLabel',{'Success', 'False Success', 'Failure'})
set(gca,'YLim',[0,1]);

SindZ = find(SuccessZ);

%saveas(gcf,'BarPlotZ.jpg')

%% CDF Plot

figure(3)

hold on

cdf = cdfplot(time(Sind));
set(cdf,'color','r')

cdfZ = cdfplot(timeZ(SindZ));
set(cdfZ,'color','b')

title('Cumulative Distribution Function')
legend('Gradient','Zeta','Location','SouthEast')
set(get(gca,'XLabel'),'String','Time')
set(get(gca,'YLabel'),'String','Percentage of Runs completed')

saveas(gcf,'CDFPlot.jpg')

hold off

%% Results Matrices

header = ['e_0, a_0, t_0, p1_0, p2_0, EXITFLAG, time, e, a, t, p1, p2, P1, P2, Success, False Success, Failure'];

results = zeros(n,17);
results(:,1)     = e0;
results(:,2)     = a0;
results(:,3)     = t0;
results(:,4:5)   = p0;
results(:,6)     = EXITFLAG;
results(:,7)     = time;
results(:,8)     = es(:,1);
results(:,9)     = as(:,1);
results(:,10)    = ts(:,1);
results(:,11:12) = ps;
results(:,13:14) = Ps;
results(1,15:17) = outcome;

%dlmwrite('results.csv', results, 'delimiter', ',', 'precision', 9);

outid = fopen('results.csv', 'w+');

fprintf(outid, '%s\n', header);

for i = 1:n
    if i == 1
        outLine = regexprep(num2str(results(i,1:17),9), '  *', ',');
    else
        outLine = regexprep(num2str(results(i,1:14),9), '  *', ',');
    end
    
    fprintf(outid, '%s\n', outLine);
end

fclose(outid);


resultsZ = zeros(n,17);
resultsZ(:,1)     = e0;
resultsZ(:,2)     = a0;
resultsZ(:,3)     = t0;
resultsZ(:,4:5)   = p0;
resultsZ(:,6)     = EXITFLAGZ;
resultsZ(:,7)     = timeZ;
resultsZ(:,8)     = esZ(:,1);
resultsZ(:,9)     = asZ(:,1);
resultsZ(:,10)    = tsZ(:,1);
resultsZ(:,11:12) = psZ;
resultsZ(:,13:14) = PsZ;
resultsZ(1,15:17) = outcomeZ;

%dlmwrite('resultsZ.csv', resultsZ, 'delimiter', ',', 'precision', 9);

outid = fopen('resultsZ.csv', 'w+');

fprintf(outid, '%s\n', header);

for i = 1:n
    if i == 1
        outLine = regexprep(num2str(resultsZ(i,1:17),9), '  *', ',');
    else
        outLine = regexprep(num2str(resultsZ(i,1:14),9), '  *', ',');
    end
    
    fprintf(outid, '%s\n', outLine);
end

fclose(outid);
