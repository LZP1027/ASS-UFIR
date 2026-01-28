function [yaw_est] = yaw_calibration(yaw_est, yaw_ture)

n = size(yaw_est, 2);

% % for j = 1:n
% %     if yaw_est(j) < 0
% %         yaw_est(j) = yaw_est(j) + 360;
% %     else
% %         yaw_est(j) = yaw_est(j);
% %     end
error = abs(yaw_ture - yaw_est);
for i = 1:n
     if error(i) > 180
         % % % if yaw_est(i) > 0
         % % %     yaw_est(i) = yaw_est(i) - 354;
         % % % elseif yaw_est(i) < 0
         % % %     yaw_est(i) = yaw_est(i) + 354;
         % % % end
         yaw_est(i) = -yaw_est(i);
     end
end
end 