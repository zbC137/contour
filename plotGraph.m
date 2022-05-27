L1 = generate_L(15, 2);
L2 = generate_L(8, 2);

figure(1)
set(gcf, 'Position', [300, 100, 720, 312], 'Color', 'White');

subplot('Position', [0.05, 0.1, 0.425, 0.77]);
plot(graph(-L1+diag(diag(L1))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
'MarkerSize', 5, 'NodeFontSize', 12);
title('(a)');

subplot('Position', [0.525, 0.1, 0.425, 0.77]);
plot(graph(-L2+diag(diag(L2))), 'LineWidth', 1.5,  'EdgeColor', 'r',...
'MarkerSize', 5, 'NodeFontSize', 12);
title('(b)');
