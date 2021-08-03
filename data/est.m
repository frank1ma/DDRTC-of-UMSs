final_alpha1(:,1)=alpha_dig_case1(:,1);
final_alpha1(:,2)=alpha_dig_case1(:,2)/max(alpha_dig_case1(:,2))*max(data_alpha(:,2));
figure;
plot(final_alpha1(:,1),final_alpha1(:,2))

final_theta1(:,1)=theta_dig_case1(:,1);
final_theta1(:,2)=theta_dig_case1(:,3)/max(theta_dig_case1(:,3))*max(data_theta(:,3));
figure;
plot(final_theta1(:,1),final_theta1(:,2))

