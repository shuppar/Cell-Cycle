%% Shuppar script for segregating cells in different stages of cell cycle

function out = CellCycleStage(FileName)

%% 
a = load(FileName);

if exist('G1_peak.dat', 'file')
    
    peaks = load('G1_peak.dat');
    g1 = peaks(1,1);
    s1 = peaks(1,3);
    g2 = peaks(1,2);
    s2 = peaks(1,4);
    

else
            fprintf('Lemme fix this.\n');
            bimodefit(FileName);
            peaks = load('G1_peak.dat');
            g1 = peaks(1,1);
            s1 = peaks(1,3);
            g2 = peaks(1,2);
            s2 = peaks(1,4);
    
end

for i = 1:size(a,1)
   
    if (a(i,1) < (g1-3*s1))
        
       f = fopen('subG1.dat','a');  
       fprintf(f,'%d\t', a(i,:));
       fprintf(f,'\n');
       
    else
        if (a(i,1) < (g1+2*s1))          % More than 83.5 % of the first gaussian is defined as G1 part with this.
        
            f = fopen('G1.dat','a');  
            fprintf(f,'%d\t', a(i,:));
	    fprintf(f,'\n');
             
            else
                if (a(i,1) < (g2-0.80*s2))  % More than 83.5 % of the second gaussian is defined as G2/M part with this.
        
                    f = fopen('S.dat','a');  
                    fprintf(f,'%d\t', a(i,:));
 		    fprintf(f,'\n');
                    
                else
                    if (a(i,1) < (g2+3*s2))
        
                        f = fopen('G2_M.dat','a');  
                        fprintf(f,'%d\t', a(i,:));
			fprintf(f,'\n');
       
                    else
                                                    
                            f = fopen('BeyondG2_M.dat','a');  
                            fprintf(f,'%d\t', a(i,:));
			    fprintf(f,'\n');
        
                    end
                end
        end
        
    end
    
    
end


%% For plotting bar graph with Control and aphidicolin combined. Not to run with this code, and if you want to, you have to make necessary changes.

% Now the point is after you have fed in the values in a, b, c, d and as, bs, cs, ds. Just copy the code hereafter and paste it in your command window with necessary changes.


% a = mean(load('G1.dat'));
% as = std(load('G1.dat')); % and so on...

% x = [a(2) b(2) c(2); d(2) e(2) f(2)];             % a is control G1.
% xs = [as(2) bs(2) cs(2); ds(2) es(2) fs(2)];      % d is Aphidic G1.
% y = [a(3) b(3) c(3); d(3) e(3) f(3)];             % b is Control S.
% ys = [as(3) bs(3) cs(3); ds(3) es(3) fs(3)];      % e is Aphidic S, and C is Control G2/M and f is Aphidic G2/m.
% 



% bb=bar(x'); hold all
% for i = 1:3  
%     j = 1:2;                           % i = number of groups and j = number of observations in each group.
%     t = -0.35 + i + 1/4 * j;           % You might have to tweak the values so as to get the errorbars at proper position.
%     errorbar(t, x(j,i), xs(j,i), '.'); 
% end


%% For plotting box plots
 
% 1) in case you want to use the inbuilt matlab function:
% 
% v = load('G1.dat');
% W = load('S.dat');
% x = load('G2_M.dat');
% a = v(:,3);              % Cyto to nuc ratio, p53 for G1 cells
% b = w(:,3);
% c = x(:,3);
% d = [a, b, c];
% e = ones([1, length(a)]);
% f = 2*ones([1, length(b)]);
% g = 3*ones([1, length(c)]);
% h = [e, f, g];
% clear a b c e f g;
% boxplot(d, h);
% 
% 
% 2) or else if you want to use some user defined thing like bplot:
% 
% v = load('G1.dat');
% W = load('S.dat');
% x = load('G2_M.dat');
% a = v(:,3);              % Cyto to nuc ratio, p53 for G1 cells
% b = w(:,3);
% c = x(:,3);
% bplot(a);
% hold on
% bplot(b, 3);
% bplot(c, 5);
% hold off;

%%

