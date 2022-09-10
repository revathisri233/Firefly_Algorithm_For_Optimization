%% Step2: Defining Objective Function

function out = fun(x)
x1=x(:,1);
x2=x(:,2);
out = x1.^2+x2.^2;
end