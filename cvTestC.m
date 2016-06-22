clc;clearvars;
%test para init
iter=1; % number of tests
[da_out, da_out1, da_out2,...
 err, err1, err2]=deal([]); % container for monitored outputs

s = 6;
abstol =1e-3;
%hyperbox = [-Inf*ones(1,s) ; sqrt(s)*rand(1,s)];
%l = -sqrt(s)*rand(1,s);u = sqrt(s)*rand(1,s);
l = -2*ones(1,s);u = ones(1,s);
hyperbox = [l;u];
sig = 0.6;  mu=0; cov = sig*ones(s,s);cov(1:s+1:s*s) = 1;

f = @(x) exp(-0.5*sum(x.*(cov\x')',2))/(sqrt(prod(eig(cov))*(2*pi)^s));
exactf = mvncdf(l,u,mu,cov);

%{
a = (1-exp(-l.^2/2))./(1-l.^2); b= 1-a;
poly = @(x) prod(repmat(a,size(x,1),1).*x.^2 + repmat(b,size(x,1),1),2); 
polymu = prod( (1/3*a.*u.^3 + b.*u)-(1/3*a.*l.^3 + b.*l) ); 
%}
fcv = @(x) exp(-sum(abs(x),2)); 
fcvmu = prod(1-exp(-u)+1-exp(l)); 
f1.func = {f,fcv};
f1.cv = fcvmu;

% begin testing 
for i=1:iter
    [q,out] = cubSobolcv_g(f,hyperbox,'uniform',abstol,0);
    err = [err;abs(q-exactf)];
    da_out = [da_out;[out.time, out.n, out.exitflag]];
end
fprintf('\n Results of cubSobolcv_g: \n');
fprintf('q=%.10f  \n', q);
%fprintf('exact=%.10f  \n', exactf);
fprintf('ave err(|q-exact|): %.10f  \n', mean(err));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out(:,2)));
fprintf('accumulated exitflag: %4d \n', sum(da_out(:,3)));

% begin testing with single variate 
for i=1:iter
    [q1,out1] = cubSobolcv_g(f1,hyperbox,'uniform',abstol,0);
    err1 = [err1;abs(q1-exactf)];
    da_out1 = [da_out1;[out1.time, out1.n, out1.exitflag]];
end
fprintf('\n Results of cubSobolcv_g(single cv): \n');
fprintf('q1=%.10f  \n',q1);
fprintf('avg err(|q1-exact|): %.10f  \n', mean(err1));
fprintf('avg time of cubSobolcv_g: %.4f \n', mean(da_out1(:,1))); 
fprintf('avg n of cubSobolcv_g: %.2f \n', mean(da_out1(:,2)));
fprintf('accumulated exitflag: %4d \n', sum(da_out1(:,3)));
fprintf('beta: %.4f \n', out1.beta);




