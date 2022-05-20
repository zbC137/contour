function plotGif(flag, data, N)

pic_num = 1;

for i = 1:50:size(data.t, 2)
% real time evolution
figure(1);
set(gcf,'Position', [200, 200, 1000, 300], 'Color', 'White');

subplot('Position', [0.0475, 0.2, 0.27, 0.7])
plot(data.f1.cd(1:2:end, i), data.f1.cd(2:2:end, i), 'k',...
        data.f1.cds(1:2:end, i), data.f1.cds(2:2:end, i), 'r',...
        data.f1.cs(1:2:end, i), data.f1.cs(2:2:end ,i),'g',...
        data.f1.c(1:2:end,  i), data.f1.c(2:2:end,  i), 'b*', 'LineWidth', 1);
% grid on;
axis([-10, 25,-12, 23]);
title(['t = ', num2str(data.t(i)), 's']);
% legend('Original Curve', 'Re-paramerterized Curve', 'Real-time Curve',...
%     'Agents Real-time Positions',  'Location', 'northwest');
% legend('boxoff')
% drawnow;

% F1 = getframe(gcf);
% I1 = frame2im(F1);
% [I1, map1] = rgb2ind(I1, 256);
% 
% if pic_num == 1
%     imwrite(I1, map1, ['pic\new contour\', flag, '.gif'], 'gif', ...
%         'LoopCount', inf, 'DelayTime', 0);
% else
%     imwrite(I1, map1, ['pic\new contour\', flag, '.gif'], 'gif', ...
%         'WriteMode', 'append', 'DelayTime', 0);
% end

% real time coff err
% figure(12)
subplot('Position',  [0.365, 0.2, 0.27, 0.7])
for j = 1:4*N+2
    plot(data.t(1:i), data.dcoff(j, 1:i), 'LineWidth', 1, 'Color', 'k'); hold on;
end
% grid on;
xlabel('time(s)'); ylabel('coefficient errors');
title(['t = ', num2str(data.t(i)), 's']);
% drawnow;

% F2 = getframe(gcf);
% I2 = frame2im(F2);
% [I2, map2] = rgb2ind(I2, 256);
% 
% if pic_num == 1
%     imwrite(I2, map2, ['pic\new contour\', flag, '_coff_err.gif'], 'gif', ...
%         'LoopCount', inf, 'DelayTime', 0);
% else
%     imwrite(I2, map2, ['pic\new contour\', flag, '_coff_err.gif'], 'gif', ...
%         'WriteMode', 'append', 'DelayTime', 0);
% end

% real time pos err
% figure(13)
subplot('Position',  [0.6825, 0.2, 0.27, 0.7])
plot(data.t(1:i), data.e(:, 1:i), 'LineWidth', 1, 'Color', 'k');
% grid on;
xlabel('time(s)'); ylabel('position errors');
title(['t = ', num2str(data.t(i)), 's']);

drawnow;

F3 = getframe(gcf);
I3 = frame2im(F3);
[I3, map3] = rgb2ind(I3, 256);

if pic_num == 1
    imwrite(I3, map3, ['pic\new contour\', flag, '.gif'], 'gif', ...
        'LoopCount', inf);
else
    imwrite(I3, map3, ['pic\new contour\', flag, '.gif'], 'gif', ...
        'WriteMode', 'append');
end

pic_num = pic_num + 1;

end