function mask = CellCycleMask(DAPI, L, peak1)
%%

g1 = peak1(1,1);
s1 = peak1(1,3);
g2 = peak1(1,2);
s2 = peak1(1,4);
Red = zeros(size(L));
Blue = Red; Green = Red;
NumLab = regionprops(L, DAPI, 'PixelIdxList', 'MeanIntensity', 'Area');
[~, num] = bwlabel(L);

%%

for i = 1:num
%         
        if (NumLab(i).MeanIntensity*NumLab(i).Area < (g1-4*s1))  
            Red(NumLab(i).PixelIdxList) = .6;                            % Yellow for sub-G1.
            Green(NumLab(i).PixelIdxList) = .8;
            Blue(NumLab(i).PixelIdxList) = .5;
            fprintf('G1-\n');
          
        else
             if (NumLab(i).MeanIntensity*NumLab(i).Area < (g1+2.0*s1))
                 Red(NumLab(i).PixelIdxList) = .6;                       % Red for G1
                 Green(NumLab(i).PixelIdxList) = .8;
                 Blue(NumLab(i).PixelIdxList) = .5;
                 fprintf('G1\n');
                 
             else
                    if (NumLab(i).MeanIntensity*NumLab(i).Area < (g2-0.8*s2))
                        Red(NumLab(i).PixelIdxList) = .3;                  % Green for S.
                        Green(NumLab(i).PixelIdxList) = .5;
                        Blue(NumLab(i).PixelIdxList) = .9;
                        fprintf('S\n');
                    
                    else
                        if (NumLab(i).MeanIntensity*NumLab(i).Area < (g2+3*s2))
                            Red(NumLab(i).PixelIdxList) = .9;              % Blue for G2/M.
                            Green(NumLab(i).PixelIdxList) = .7;
                            Blue(NumLab(i).PixelIdxList) = .3;  
                            fprintf('G2\n');
                            
                        else
                            
                                 Red(NumLab(i).PixelIdxList) = .9;         % Green-blue for beyond M.
                                 Green(NumLab(i).PixelIdxList) = .7;
                                 Blue(NumLab(i).PixelIdxList) = .3;
                                 fprintf('G2+');
                        
                        end
                        
                    end
                    
              end
             
             
        end
        
        
end

mask = cat(3, Red, Green, Blue);
figure, imshow(mask, []);


%%




