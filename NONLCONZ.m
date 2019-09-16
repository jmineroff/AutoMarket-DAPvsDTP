function [C,Ceq,DC,DCeqt] = NONLCONZ(x,e,a,t,w,f,c,ACoef,UCoef,cCoef)

C = [];

e(1) = x(1);
a(1) = x(2);
t(1) = x(3);
p(1) = x(4);
p(2) = x(5);

% Vehicle Attributes
Ceq(1) = 1000 ./ ( e(1) - ACoef(1) ) ...
    - ( ACoef(2) ...
    + ACoef(3) .* exp( -a(1) ) ...
    + ACoef(4) .* w(1)...
    + ACoef(5) .* w(1).* a(1) ...
    + ACoef(6) .* t(1) ...
    + ACoef(7) .* t(1) .* a(1)^2 );

U(1) = UCoef(1) * p(1)...   % 10,000 2006$
    + UCoef(2) / e(1)... % mpg
    + UCoef(3) / a(1)... % s
    + UCoef(4) * f(1)... % 10,000 in^2
    + UCoef(5);

U(2) = UCoef(1) * p(2)...   % 10,000 2006$
    + UCoef(2) / e(2)... % mpg
    + UCoef(3) / a(2)... % s
    + UCoef(4) * f(2)... % 10,000 in^2
    + UCoef(5);

P(1) = exp( U(1) ) / ( 1.0 + exp( U(1) ) + exp( U(2) ) );
P(2) = exp( U(2) ) / ( 1.0 + exp( U(1) ) + exp( U(2) ) );

c(1) = cCoef(1)...
    + cCoef(2) * exp( -a(1) )...
    + cCoef(3) * t(1) ...
    + cCoef(4) * w(1) ...
    + cCoef(5) * w(1) * a(1);

DpU = UCoef(1);

m(1) = p(1) - c(1);
m(2) = p(2) - c(2);

pi(1) = P(1)*m(1);
pi(2) = P(2)*m(2);

z(1) = pi(1) - 1/DpU;
z(2) = pi(2) - 1/DpU;

phi(1) = m(1) - z(1);
phi(2) = m(2) - z(2);

% Price Equilibrium 1 Zeta
Ceq(2) = phi(1);

% Price Equilibrium 2 Zeta
Ceq(3) = phi(2);

if nargout > 2
    
    DeU = -UCoef(2)/( e(1)^2 );
    DaU = -UCoef(3)/( a(1)^2 );
    DtU = 0;
    
    Dec = 0;
    Dac = -cCoef(2)*exp( -a(1) ) + cCoef(5)*w(1);
    Dtc = cCoef(3);
    
    omP = 1 - P(1);
    
    DC = [];
    
    % Jacobian
    DCeq = zeros(3,5);
    
    % Design constraint
    DCeq(1,1) = -1000/((e(1) - ACoef(1))^2);
    DCeq(1,2) = ACoef(3)*exp(-a(1)) ...
        - ACoef(5)*w(1) ...
        - 2*ACoef(7)*a(1)*t(1);
    DCeq(1,3) = -ACoef(6) - ACoef(7)*a(1)^2;
    
    % Price equilibrium 1 Zeta
    DCeq(2,1) = -DeU*P(1)*omP *m(1) + ( P(1) - 1 )*Dec;
    DCeq(2,2) = -DaU*P(1)*omP *m(1) + ( P(1) - 1 )*Dac;
    DCeq(2,3) = -DtU*P(1)*omP *m(1) + ( P(1) - 1 )*Dtc;
    DCeq(2,4) = -DpU*P(1)*omP *m(1) - ( P(1) - 1 );
    DCeq(2,5) =  DpU*P(1)*P(2)*m(1);
    
    % Price equilibrium 2 Zeta
    DCeq(3,1) = DeU*P(1)*P(2)*m(2);
    DCeq(3,2) = DaU*P(1)*P(2)*m(2);
    DCeq(3,3) = DtU*P(1)*P(2)*m(2);
    DCeq(3,4) = DpU*P(1)*P(2)*m(2);
    DCeq(3,5) = 1 - DpU*P(2)*phi(2);
    
    DCeqt = DCeq';
end

% global Safety;
% 
% if Safety == 1
%     
%     Safety = 0;
%     
%     CheckGrad(@(x)NONLCONZ(x,e,a,t,w,f,c,ACoef,UCoef,cCoef), x, 3);
%     
%     Safety = 1;
%     
% end

end
