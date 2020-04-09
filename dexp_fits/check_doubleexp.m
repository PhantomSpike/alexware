%check gradient
Ainit =  (0.5+rand)./100;
Binit = (0.5+rand)./1000;
alphainit = rand;
betainit = 4+4.*rand;
%Fake data
T = randn(25,1)./100000 + [1 -0.1 -0.08 -0.06 -0.04 -0.02 -0.01 -0.01 -0.005 -0.005 zeros(1,15)]'./1000;

vinit = [Ainit; Binit; alphainit; betainit];
L = (0:1:24)';


checkgrad('doubleexp_E_dEdv',vinit,10e-5,L,T);

%fit
for n = 1:100;
    Ainit =  (0.5+rand)./1000;
    Binit = (0.5+rand)./1000;
    alphainit = rand;
    betainit = 10.*rand;
    
    vinit = [Ainit; Binit; alphainit; betainit];
    
    [vopt(:,n),Fopt(n)] = minFunc(@doubleexp_E_dEdv,vinit,[],L,T);
end

[~,best] = min(Fopt);
vopt = vopt(:,best);

plot(L,T,'k.');

hold on

A = vopt(1);
B = vopt(2);
alpha = vopt(3);
beta = vopt(4);

expLal = exp(-L./abs(alpha));
expLbe = exp(-L./abs(beta));
F = abs(A).*expLal - abs(B).*expLbe;

plot(L,F,'r')