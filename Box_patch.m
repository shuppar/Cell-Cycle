


%%
cd Control/Data
% CellCycleStage('Intensity.dat');
A4 = load('Intensity.dat');
A1 = load('G1.dat'); 
A2 = load('S.dat'); 
A3 = load('G2_M.dat');
A11m = mean(A1(:,3)); A21m = mean(A2(:,3)); A31m = mean(A3(:,3)); A41m = mean(A4(:,3));
A12m = mean(A1(:,8)); A22m = mean(A2(:,8)); A32m = mean(A3(:,8)); A42m = mean(A4(:,8));
% A1e = std(A1(:,3))/sqrt(length(A1)); A2e = std(A2(:,3))/sqrt(length(A2)); A3e = std(A3(:,3))/sqrt(length(A3)); A4e = std(A4(:,3))/sqrt(length(A4));

%%
cd ../..
cd 1/Data
% CellCycleStage('Intensity.dat');
B1 = load('G1.dat');
B2 = load('S.dat');
B3 = load('G2_M.dat');
B4 = load('Intensity.dat');
B11m = mean(B1(:,3)); B21m = mean(B2(:,3)); B31m = mean(B3(:,3)); B41m = mean(B4(:,3));
B12m = mean(B1(:,8)); B22m = mean(B2(:,8)); B32m = mean(B3(:,8)); B42m = mean(B4(:,8));
% B1e = std(B1(:,3))/sqrt(length(B1)); B2e = std(B2(:,3))/sqrt(length(B2)); B3e = std(B3(:,3))/sqrt(length(B3)); B4e = std(B4(:,3))/sqrt(length(B4));

%%
cd ../..
cd 3/Data
% CellCycleStage('Intensity.dat');
C1 = load('G1.dat');
C2 = load('S.dat');
C3 = load('G2_M.dat');
C4 = load('Intensity.dat');
C11m = mean(C1(:,3)); C21m = mean(C2(:,3)); C31m = mean(C3(:,3)); C41m = mean(C4(:,3));
C12m = mean(C1(:,8)); C22m = mean(C2(:,8)); C32m = mean(C3(:,8)); C42m = mean(C4(:,8));
% C1e = std(C1(:,3))/sqrt(length(C1)); C2e = std(C2(:,3))/sqrt(length(C2)); C3e = std(C3(:,3))/sqrt(length(B3)); C4e = std(C4(:,3))/sqrt(length(C4));

%%
cd ../..
cd 10/Data
% CellCycleStage('Intensity.dat');
D1 = load('G1.dat');
D2 = load('S.dat');
D3 = load('G2_M.dat');
D4 = load('Intensity.dat');
D11m = mean(D1(:,3)); D21m = mean(D2(:,3)); D31m = mean(D3(:,3)); D41m = mean(D4(:,3));
D12m = mean(D1(:,8)); D22m = mean(D2(:,8)); D32m = mean(D3(:,8)); D42m = mean(D4(:,8));
% D1e = std(D1(:,3))/sqrt(length(D1)); D2e = std(D2(:,3))/sqrt(length(D2)); D3e = std(D3(:,3))/sqrt(length(D3)); D4e = std(D4(:,3))/sqrt(length(D4));

%%
cd ../..
cd 30/Data
% CellCycleStage('Intensity.dat');
E1 = load('G1.dat');
E2 = load('S.dat');
E3 = load('G2_M.dat');
E4 = load('Intensity.dat');
E11m = mean(E1(:,3)); E21m = mean(E2(:,3)); E31m = mean(E3(:,3)); E41m = mean(E4(:,3));
E12m = mean(E1(:,8)); E22m = mean(E2(:,8)); E32m = mean(E3(:,8)); E42m = mean(E4(:,8));
% E1e = std(E1(:,3))/sqrt(length(E1)); E2e = std(E2(:,3))/sqrt(length(E2)); E3e = std(E3(:,3))/sqrt(length(E3)); E4e = std(E4(:,3))/sqrt(length(E4));

%% Normalizing with respect to the maximum of means of respective columns of all ...
% the phases for all the doses

N1 = max([A11m, B11m, C11m, D11m, E11m,A21m, B21m, C21m, D21m, E21m,A31m, B31m, C31m, D31m, E31m]);
N2 = max([A12m, B12m, C12m, D12m, E12m,A22m, B22m, C22m, D22m, E22m,A32m, B32m, C32m, D32m, E32m]);


%% For G1 phase cells

X1 = [A11m, B11m, C11m, D11m, E11m]/N1; Y1 = [A12m, B12m, C12m, D12m, E12m]/N2;
Nt1 = X1(1)+Y1(1);
X1 = 10*X1/Nt1; Y1 = 10*Y1/Nt1;
Y1 = 10-Y1;

% v11 = [1 10; 1 Y1(1); 2 Y1(1); 2 Y1(2); ...
%        3 Y1(2); 3 Y1(3); 4 Y1(3); 4 Y1(5); 5 Y1(5); 5 10; ...
%        1 0; 1 X1(1); 2 X1(1); 2 X1(2); ...
%        3 X1(2); 3 X1(3); 4 X1(3); 4 X1(5); 5 X1(5); 5 0];
% f11 = [1 2 3 4 5 6 7 8 9 10; 11 12 13 14 15 16 17 18 19 20];


v11 = [1 10; 1 Y1(1); 2 Y1(1); 2 Y1(2); ...
       3 Y1(2); 3 Y1(3); 4 Y1(3); 4 Y1(5); 5 Y1(5); 5 10];
f11 = [1 2 3 4 5 6 7 8 9 10];
figure, p11 =  patch('Faces', f11, 'Vertices', v11, 'FaceColor', [254,246,0]/255,  ...
    'EdgeColor', [254,216,0]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p11.LineStyle = 'none';

v12 = [1 0; 1 X1(1); 2 X1(1); 2 X1(2); ...
       3 X1(2); 3 X1(3); 4 X1(3); 4 X1(5); 5 X1(5); 5 0];
f12 = [1 2 3 4 5 6 7 8 9 10];
p12 =  patch('Faces', f12, 'Vertices', v12, 'FaceColor', [254,246,0]/255, ...
    'EdgeColor', [254,216,0]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p12.LineStyle = 'none';
hatchfill2(p11,'FaceColor',[254,246,0]/255,'HatchStyle','single','HatchAngle',60,'HatchDensity',13, 'HatchColor', 'black');
hatchfill2(p12,'FaceColor',[254,246,0]/255,'HatchStyle','single','HatchAngle',120,'HatchDensity',13, 'HatchColor', 'black');

xlim([-5 10])
ylim([-5 15])

hold on
plot([1, 1, 5, 5, 1], [0, 10, 10, 0, 0], 'LineWidth', 3, 'Color', [254,216,0]/255)
hold off

%% For S phase cells

X2 = [A21m, B21m, C21m, D21m, E21m]/N1; Y2 = [A22m, B22m, C22m, D22m, E22m]/N2;
Nt2 = X2(1)+Y2(1);
X2 = 10*X2/Nt2; Y2 = 10*Y2/Nt2;
Y2 = 10-Y2;

v21 = [1 10; 1 Y2(1); 2 Y2(1); 2 Y2(2); ...
       3 Y2(2); 3 Y2(3); 4 Y2(3); 4 Y2(5); 5 Y2(5); 5 10];
f21 = [1 2 3 4 5 6 7 8 9 10];
figure, p21 =  patch('Faces', f21, 'Vertices', v21, 'FaceColor', [250,155,0]/255,  ...
    'EdgeColor', [235,110,0]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p21.LineStyle = 'none';

v22 = [1 0; 1 X2(1); 2 X2(1); 2 X2(2); ...
       3 X2(2); 3 X2(3); 4 X2(3); 4 X2(5); 5 X2(5); 5 0];
f22 = [1 2 3 4 5 6 7 8 9 10];
p22 =  patch('Faces', f22, 'Vertices', v22, 'FaceColor', [250,155,0]/255, ...
    'EdgeColor', [235,110,0]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p22.LineStyle = 'none';
hatchfill2(p21,'FaceColor',[250,155,0]/255,'HatchStyle','single','HatchAngle',60,'HatchDensity',13, 'HatchColor', 'black');
hatchfill2(p22,'FaceColor',[250,155,0]/255,'HatchStyle','single','HatchAngle',120,'HatchDensity',13, 'HatchColor', 'black');

xlim([-5 10])
ylim([-5 15])

hold on
plot([1, 1, 5, 5, 1], [0, 10, 10, 0, 0], 'LineWidth', 3, 'Color', [235,110,0]/255)
hold off


%% For G2/M phase cells

X3 = [A31m, B31m, C31m, D31m, E31m]/N1; Y3 = [A32m, B32m, C32m, D32m, E32m]/N2;
Nt3 = X3(1)+Y3(1);
X3 = 10*X3/Nt3; Y3 = 10*Y3/Nt3;
Y3 = 10-Y3;

v31 = [1 10; 1 Y3(1); 2 Y3(1); 2 Y3(2); ...
       3 Y3(2); 3 Y3(3); 4 Y3(3); 4 Y3(5); 5 Y3(5); 5 10];
f31 = [1 2 3 4 5 6 7 8 9 10];
figure, p31 =  patch('Faces', f31, 'Vertices', v31, 'FaceColor', [150,40,20]/255,  ...
    'EdgeColor', [80,17,6]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p31.LineStyle = 'none';

v32 = [1 0; 1 X3(1); 2 X3(1); 2 X3(2); ...
       3 X3(2); 3 X3(3); 4 X3(3); 4 X3(5); 5 X3(5); 5 0];
f32 = [1 2 3 4 5 6 7 8 9 10];
p32 =  patch('Faces', f32, 'Vertices', v32, 'FaceColor', [150,40,20]/255, ...
    'EdgeColor', [80,17,6]/255, 'LineWidth', 3, 'FaceAlpha', 0.5);
% p32.LineStyle = 'none';
hatchfill2(p31,'FaceColor',[150,40,20]/255,'HatchStyle','single','HatchAngle',60,'HatchDensity',13, 'HatchColor', 'black');
hatchfill2(p32,'FaceColor',[150,40,20]/255,'HatchStyle','single','HatchAngle',120,'HatchDensity',13, 'HatchColor', 'black');

xlim([-5 10])
ylim([-5 15])

hold on
plot([1, 1, 5, 5, 1], [0, 10, 10, 0, 0], 'LineWidth', 3, 'Color', [80,17,6]/255)
hold off

%%

% %%
% 
% Data = [A1m B1m C1m D1m E1m ; A2m B2m C2m D2m E2m; A3m B3m C3m D3m E3m]/A4m;
% Error = [A1e B1e C1e D1e E1e; A2e B2e C2e D2e E2e; A3e B3e C3e D3e E3e]/A4m;
% axes1 = axes('Parent',figure);
% CT=cbrewer('seq', 'YlOrBr', 5);
% 
% h = barwitherr(Error, Data, 'LineWidth', 2);
% set(gca,'XTick',[1 2 3],'XTickLabel',{'G1','S','G2/M'});
% fh = gcf; colormap(CT);
% 
% lbl = {'Control', '0.4 \mug/ml', '0.8 \mug/ml', '1.6 \mug/ml', '3.2 \mug/ml'};
% 
% legendflex(h, lbl, 'anchor', {'nw','nw'}, ...
% 'FontName','Times', ...
% 'FontWeight','bold', ...
% 'buffer', [5 -5], ...
% 'ncol', 3, ...
% 'fontsize', 37, ...
% 'xscale', 0.8, ...
% 'box', 'off');
% 
% 
% set(axes1,'FontName','Times','FontSize',37,'FontWeight','bold')
% title('Puro Only', 'FontSize', 37, 'FontWeight', 'bold', 'FontName', 'Times');

%% The End! :)