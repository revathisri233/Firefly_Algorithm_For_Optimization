
 %% Firefly Algorithm MATLAB code
format short
clear all
clc
%% step1: Input Parameters
D = 2;              % Dimension of the problem
lb = [-5.12 -5.12]; % Lower bound of the variables
ub = [5.12 5.12];   % Upper bound of the variables
N = 20;             % Population size
alpha =1.0;         % Randomness strength 0-1(highly random)
beta0 =1.0;         % Attractiveness constant
gamma =0.01;        % Absorption Coefficient
delta =0.97;        % Randomness Reduction Factor delta = 10^(-5/iter_max)
iter_max=100;       % Maximum number of iterations
%% Step 2: Defining objective function - Benchbark problem De-Jong's function f(x)= sum(x_i)^2 for i=1:D where D= dimensions
%function out = fun(x)
%for i=1:D
    %x(i)=x(:,i);
    %out = sum(x(i)).^2;
%end
%end
%% step3: Generate initial population randomly
for i=1:N
    for j=1:D
        pop(i,j)=lb(:,j)+rand.*(ub(:,j)-lb(:,j)); % use x=L+rand.*(U-L)
    end
end
%% Evaluate Objective Function
fx = fun(pop);     
alpha= alpha*delta;  % Reduce alpha by a factor theta
scale= abs(ub-lb);   % Scale of the problem

%% step 4:FIREFLY ALGORITHM MAIN LOOP STARTS 
for iter = 1:iter_max
    %% for n fireflies,there are a total of two loops 
    for i=1:N
        for j=1:N
            fx(i)=fun(pop(i,:));
            %% Brighter / more attractive firefly
            if fx(i) < fx(j)
                pop(i,:) = pop(i,:);
            elseif fx(i) > fx(j)
                xi=pop(i,:);
                xj=pop(j,:);
                r=sqrt(sum((xi-xj).^2));
                beta=beta0*exp(-gamma*r.^2);
                steps=alpha.*(rand(1,D)-0.5).*scale;
                xnew=xi+beta*(xj-xi)+steps;  % New position update equation
        %% step 5:Check the bounds
        for i=1:size(xnew,2)
            if xnew(i)>ub(i)
                xnew(i)=ub(i);
            elseif xnew(i)<lb(i)
                xnew(i)=lb(i);
            end
        end
        %% step6: perform greedy selection
        fnew = fun(xnew);   % for minimisation, fnew<fold,then update the solution,otherwise not
        if fnew<fx(i)
            fx(i)=fnew;
            pop(i,:)=xnew;
        end
            end % end for "if fx(i)<fx(j)"
        end % end for j
    end %end for i
    %% Memorize the best
    [optval, optind] = min(fx);
    Bestfx(iter) = optval;
    Bestx(iter,:)=pop(optind,:);
    %% Show iteration information
    disp(['Iteration' num2str(iter)...
        ':Best Cost =' num2str(Bestfx(iter))]);
    
    %% Plotting the result
    plot(Bestfx, 'LineWidth', 2);
    xlabel('Iteration number');
    ylabel('Fitness value')
    title('Convergence vs Iteration')
    grid on
end %%%% end of the ITERATION LOOP
     