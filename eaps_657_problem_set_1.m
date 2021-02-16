%PROBLEM A
x = -20:1:20;
gzobs = 41.93*(0.5)*(0.8)*(1./(x.^2./5.^2+1))
plot(x,gzobs)
title('Surface Anomaly - Problem A')
ylabel('gzobs')
xlabel('x')
axis equal

%PROBLEM B
%z should vary from 3 to 20, dp should vary from 0.1 to 2.
%increment so that z and delta_rho are 1 x 41 vectors
increment_delta_rho = (2-0.1)/40;
delta_rho = 0.1:increment_delta_rho:2;
increment_z = (20-3) / 40;
z = 3:increment_z:20;

%run for loop to iterate through each delta_rho / z combination on grid.
%For problem b, there are 41 x 41 points to iterate through
num_zeros = size(x); 
error = zeros(num_zeros(2)); %for now, create empty 41 x 41 matrix of 0s to be filled in as the for loop below progresses
for i = 1:num_zeros(2) %z 
    for j = 1:num_zeros(2) %delta_rho
        %delta_rho(j);
        %z(i);
        %error(i,j)
        % comments to help ensure code is computing as desired: fprintf('place in matrix is(%.0f , %.0f) where delta_rho and z are %f and %f\n', i, j, delta_rho(j), z(i))
        gzpre = (41.93*delta_rho(j) * (4 / z(i))) * (1./((x.^2 ./ z(i)^2) + 1)); %gzpre varies with each i and delta_rho. Each 
        % point on a grid results in a 1 x 41 vector, because x is a 1 x 41
        % vector
        sq_residual = (gzpre - gzobs).^2; %create squared residual
        error(i,j) = sqrt(sum(sq_residual)); %square root the sum of the residuals, assign that value to be the error at the corresponding i,j on the grid
        % error(i,j); %output to ensure code generates desired output. 
        gzpre = 0; %reset gzpre for next iteration of calculation - unclear if value from this iteration of for loop
        %carries over to next iteration. Just a safety net in code. 
    end
end

%find parameters with smallest error
minimum = min(min(error));
[b, a] = find(error == minimum);
best_parameters = [z(1) + (a - 1) * increment_z, delta_rho(1) + (b-1) * increment_delta_rho]

%plot 3d error surface
figure;
surf(z, delta_rho, error)
title('3d Error Surface - Problem B')
xlabel('z')
ylabel('delta rho')
zlabel('Error')

%PROBLEM C
%re-assign x to a more constrained length. Go from -20:1:20 to -10:1:10 -
%this will reduce the computation by about half in the most expensive part
%of the code. 

x = -10:1:10;

%recompute gzobs for updated x
gzobs = 41.93*(0.5)*(0.8)*(1./(x.^2./5.^2+1));

%update vectors for z and delta_rho to match dimension of gzobs
delta_rho = linspace(0.1,2,21);
increment_delta_rho_problem_c = (2-0.1)/20;
increment_z_problem_c = (20-3) / 20;
z = linspace(3,20,21);

num_zeros_problem_c = size(x);
error_problem_c = zeros(num_zeros_problem_c(2));
%repeat for loop from problem b, but with updated x
 for i = 1:num_zeros_problem_c(2) %z 
    for j = 1:num_zeros_problem_c(2) %delta_rho
        delta_rho(j);
        z(i);
        error_problem_c(i,j);
        %fprintf('place in matrix is(%.0f , %.0f) where delta_rho and z are %f and %f\n', i, j, delta_rho(j), z(i))
        gzpre = (41.93*delta_rho(j) * (4 / z(i))) * (1./((x.^2 ./ z(i)^2) + 1));
        sq_residual = (gzpre - gzobs).^2;
        error_problem_c(i,j) = sqrt(sum(sq_residual));
        %error_problem_c(i,j)
        gzpre = 0; %reset gzpre for next iteration of calculation - unclear if value from this iteration of for loop
        %carries over to next iteration. Just a safety net in code.
    end
end

%find parameters with smallest error
minimum_problem_c = min(min(error_problem_c));
[b_c, a_c] = find(error_problem_c == minimum_problem_c);
best_parameters_problem_c = [z(1) + (a_c - 1) * increment_z_problem_c, delta_rho(1) + (b_c-1) * increment_delta_rho_problem_c]

%plot 3d error surface
figure;
surf(z, delta_rho, error_problem_c)
title('3d Error Surface - Problem C')
xlabel('z')
ylabel('delta rho')
zlabel('Error')

%PROBLEM D
%add noise / error into gzobs: 1-2 mgal randomly.
%assumes random numbers are generated for surface -10 < x < 10 (i.e 21
%elements with error in gzobs
random = 1 + (2-1).*rand(21,1) %create random vector to match gzobs dimension of 1 x 21 - errors between 1-2 mgal
gzobs_problem_d = gzobs + random'; %add error to gzobs

num_zeros_problem_d = size(x);
error_problem_d = zeros(num_zeros_problem_c(2));
%repeat for loop from problem b, but with updated x and noise added to
%gzobs
 for i = 1:num_zeros_problem_d(2) %z 
    for j = 1:num_zeros_problem_d(2) %delta_rho
        %delta_rho(j)
        %z(i)
        %error_problem_d(i,j)
        %fprintf('place in matrix is(%.0f , %.0f) where delta_rho and z are %f and %f\n', i, j, delta_rho(j), z(i))
        gzpre = (41.93*delta_rho(j) * (4 / z(i))) * (1./((x.^2 ./ z(i)^2) + 1));
        sq_residual = (gzpre - gzobs_problem_d).^2;
        error_problem_d(i,j) = sqrt(sum(sq_residual));
        %error_problem_d(i,j)
        gzpre = 0; %reset gzpre for next iteration of calculation - unclear if value from this iteration of for loop
        %carries over to next iteration. Just a safety net in code.
    end
end

%find parameters with smallest error
minimum_problem_d = min(min(error_problem_d));
[b_d, a_d] = find(error_problem_d == minimum_problem_d);
best_parameters_problem_d = [z(1) + (a_d - 1) * increment_z_problem_c, delta_rho(1) + (b_d-1) * increment_delta_rho_problem_c]

%plot 3d error surface
figure;
surf(z, delta_rho, error_problem_d)
title('3d Error Surface - Problem D')
xlabel('z')
ylabel('delta rho')
zlabel('Error')

