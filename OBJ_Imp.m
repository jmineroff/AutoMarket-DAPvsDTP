function [npr,nprg] = OBJ_Imp(x,e,a,t,w,f,cCoef,UCoef)

e(1) = x(1);
a(1) = x(2);
t(1) = x(3);
p(1) = x(4);
p(2) = x(5);

U(1) = UCoef(1)* p(1)... % 10,000 2006$
    + UCoef(2)/ e(1)... % mpg
    + UCoef(3)/ a(1)... % s
    + UCoef(4)* f(1)... % 10,000 in^2
    + UCoef(5);

U(2) = UCoef(1)* p(2)... % 10,000 2006$
    + UCoef(2)/ e(2)... % mpg
    + UCoef(3)/ a(2)... % s
    + UCoef(4)* f(2)...    % 10,000 in^2
    + UCoef(5);

P(1) = exp( U(1) )/ ( 1.0 + exp( U(1) ) + exp( U(2) ) );

P(2) = exp( U(2) )/ ( 1.0 + exp( U(1) ) + exp( U(2) ) );

c(1) = cCoef(1)...
    + cCoef(2) * ( exp( -a(1) ) ) ...
    + cCoef(3) * t(1) ...
    + cCoef(4) * w(1) ...
    + cCoef(5) * w(1) * a(1);

npr = - P(1) * ( p(1) - c(1) );

if nargout == 2
    
    DpU = UCoef(1);
    DeU = -UCoef(2)/( e(1)^2 );
    DaU = -UCoef(3)/( a(1)^2 );
    DtU = 0;
    
    Dec = 0;
    Dac = -cCoef(2)*exp( -a(1) ) + cCoef(5)*w(1);
    Dtc = cCoef(3);
    
    omP = 1 - P(1);
    m = p(1) - c(1);
    nprg = zeros(1,5);
    nprg(1,1) = -DeU*P(1)*omP *m(1) + P(1)*Dec; % e
    nprg(1,2) = -DaU*P(1)*omP *m(1) + P(1)*Dac; % a
    nprg(1,3) = -DtU*P(1)*omP *m(1) + P(1)*Dtc; % t
    nprg(1,4) = -DpU*P(1)*omP *m(1) - P(1);     % p(1)
    nprg(1,5) =  DpU*P(1)*P(2)*m(1);            % p(2)
    
end

% global Safety;
% 
% if Safety == 1
%     
%     Safety = 0;
%     
%     CheckGrad(@(x)( OBJ(x,e,a,t,w,f,cCoef,UCoef)), x);
%     
%     Safety = 1;
%
% end

end