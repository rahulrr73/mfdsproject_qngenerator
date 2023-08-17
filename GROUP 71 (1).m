# Q111
% Define the size of the matrix and the rank of the approximation
n = 500;
rank_approx = 100;

% Generate a random matrix with rank r
A = randn(n, rank_approx)*randn(rank_approx, n);

% Compute the SVD of A
[U, S, V] = svd(A);

% Compute the rank-r approximation of A using its SVD
A_approx = U(:,1:rank_approx)*S(1:rank_approx,1:rank_approx)*V(:,1:rank_approx)';

% Compute the relative error between A and its rank-r approximation
error = norm(A - A_approx, 'fro')/norm(A, 'fro');

% Display the relative error
fprintf('Relative error between A and its rank-%d approximation: %e\n', rank_approx, error);

% Plot the singular values of A and their cumulative sum
singular_values = diag(S);
cumulative_sum = cumsum(singular_values.^2)/sum(singular_values.^2);
figure;
subplot(2,1,1);
semilogy(singular_values, 'ko', 'LineWidth', 2);
ylabel('Singular values');
xlabel('Index');
title('Singular values of A');
subplot(2,1,2);
plot(cumulative_sum, 'k-', 'LineWidth', 2);
ylim([0,1]);
ylabel('Cumulative sum');
xlabel('Number of singular values');
title('Cumulative sum of singular values of A');

% Compute and plot the relative error as a function of the number of singular values used
max_rank = min(n, rank_approx*5);
error_rank = zeros(max_rank, 1);
for r = 1:max_rank
    A_approx = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    error_rank(r) = norm(A - A_approx, 'fro')/norm(A, 'fro');
end
figure;
semilogy(error_rank, 'k-', 'LineWidth', 2);
xlabel('Number of singular values');
ylabel('Relative error');
title('Relative error of rank-r approximations of A');

#QQ2
n = randi([50, 200]); % generate a random matrix size between 50 and 200
A = randn(n); % generate a random matrix of size nxn
b = randn(n, 1); % generate a random column vector of size nx1
tol = 1e-6; % set the tolerance for convergence
max_iters = 1000; % set the maximum number of iterations

% initialize x with a random vector of size nx1
x = randn(n, 1);

% perform the shifted power method to find the dominant eigenvalue and eigenvector of A
for k = 1:max_iters
    % compute the Rayleigh quotient of x and A to estimate the dominant eigenvalue
    lambda = (x'*A*x)/(x'*x);
    
    % subtract the estimate of the dominant eigenvalue from A
    Ap = A - lambda*eye(n);
    
    % perform one iteration of the power method on Ap
    y = Ap\x;
    
    % normalize y to have unit norm
    x = y/norm(y);
    
    % check for convergence
    if norm(Ap*x - lambda*x) < tol
        break;
    end
end

% output the dominant eigenvalue and eigenvector
fprintf('The dominant eigenvalue is %f\n', lambda);
fprintf('The corresponding eigenvector is:\n');
disp(x);
##Q5
%QQ5
%computing the probability of observing a sequence of five consecutive heads when tossing a biased coin that lands on heads with probability p. 
num_trials = 10000; % number of trials to run
num_tosses = 20; % maximum number of tosses to make
p_heads = 0.443; % probability of heads

count_success = 0; % counter for number of successful trials
success_idx = zeros(num_trials, 1); % index of successful trials

for i = 1:num_trials
    sequence = ''; % initialize sequence of tosses
    for j = 1:num_tosses
        if rand() < p_heads
            sequence = strcat(sequence, 'H'); % heads
        else
            sequence = strcat(sequence, 'T'); % tails
        end
        
        % check for three consecutive heads
        if length(sequence) >= 5 && strcmp(sequence(end-4:end), 'HHHHH')
            count_success = count_success + 1;
            success_idx(i) = j; % record index of successful trial
            break; % stop tossing if success
        end
    end
end

prob_success = count_success / num_trials % estimated probability of success

% plot the results
figure;
histogram(success_idx(success_idx>0), 'Normalization', 'probability');
xlim([0 num_tosses]);
xlabel('Number of tosses until success');
ylabel('Probability');
title('Distribution of number of tosses until three consecutive heads');

% print the estimated probability of success
fprintf('Estimated probability of success: %f\n', prob_success);


%Question6
%Write a MATLAB code to perform the following tasks:
%Plot the points on a scatter plot.
%Fit a circle to the points using the least squares method.
%Plot the circle on the same scatter plot.
%Compute the area of the circle.
%Compute the average distance between the points and the circle.
% Generate random x and y values with 1000 samples each
x = randn(3000, 1);
y = randn(3000, 1);

% Plot the points on a scatter plot
scatter(x, y, 'filled');
hold on;

% Fit a circle to the points using the least squares method
A = [x, y, ones(length(x), 1)];
b = -x.^2 - y.^2;
c = pinv(A) * b;
xc = -0.5*c(1);
yc = -0.5*c(2);
r = sqrt((xc - x(1))^2 + (yc - y(1))^2);

% Plot the circle on the same scatter plot
theta = linspace(0, 2*pi);
plot(xc + r*cos(theta), yc + r*sin(theta), 'r');

% Add labels for center and radius
text(xc + r*cos(pi/4), yc + r*sin(pi/4), sprintf('r = %.2f', r), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(xc, yc, sprintf('(%.2f, %.2f)', xc, yc), 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Compute the area of the circle
area = pi*r^2;

% Compute the average distance between the points and the circle
distances = sqrt((x - xc).^2 + (y - yc).^2);
avg_distance = mean(abs(distances - r));

% Display results
fprintf('Center: (%.2f, %.2f)\n', xc, yc);
fprintf('Radius: %.2f\n', r);
fprintf('Area: %.2f\n', area);
fprintf('Average distance: %.2f\n',avg_distance);
