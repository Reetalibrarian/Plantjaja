%%
% clearvars -except ff exec_dir upload_dir dir_list phyllotaxy full_path
clearvars -except ff exec_dir phyllotaxy full_path upload_dir

    IM = imread([full_path]);

    n = 500; % resize
    n_resize = n;
    IM0 = imresize(IM, [n NaN]);
    IMi = (double(IM0)./255).^1*255;
    maxR = max(max(IMi(:,:,1)));
    maxG = max(max(IMi(:,:,2)));
    maxB = max(max(IMi(:,:,3)));
    IMi(:,:,1) = IMi(:,:,1)/maxR*255;
    IMi(:,:,2) = IMi(:,:,2)/maxG*255;
    IMi(:,:,3) = IMi(:,:,3)/maxB*255;
    hsvIM4 = rgb2hsv(IMi);
 
% moth ROI
    IMmask = zeros(size(hsvIM4(:,:,3))); 
    IMmask(find(IMi(:,:,3)<120)) = 1;
    %figure, imshow(IMmask,[])
    %figure(444), imshow(BW,[]);
    maskofIMi = IMmask;

    BW_jieu = zeros(size(IMmask));
    for i=(size(IMmask,1)-1):-1:2
        for j=(size(IMmask,2)-1):-1:2
           if sum(sum(IMmask(i-1:i+1,j-1:j+1)))>1;
               BW_jieu(i,j)=1;
            else BW_jieu(i,j)=0;
           end
        end
    end
    
    c = 0.5*n;
    r = min(find(BW_jieu(c,:)>0));
    
    contour = bwtraceboundary(BW_jieu,[c r],'S',8,inf,'counterclockwise');
    
    %[num2str(ff) '---' dir_list(ff).name]
    
    c1Contour{ff,:} = contour(:,1);
    c2Contour{ff,:} = contour(:,end);

% smooth = 5,5
    y1 = smooth(c1Contour{ff,:},5);
    y2 = smooth(c2Contour{ff,:},2);
    [ymax1,imax1,ymin1,imin1] = extrema(y1);
    [ymax2,imax2,ymin2,imin2] = extrema(y2);
    n = size(ymax1,1)+size(ymin1,1);
    n2 = size(ymax2,1)+size(ymin2,1);
    M1_value(ff,1)=n;
    M2_value(ff,1)=n2;

% smooth = 10, 10
    smoothPar = fix(size(contour,1)*0.6);
    sy1 = smooth(c1Contour{ff,:},50);
    sy2 = smooth(c2Contour{ff,:},50);
    %sy1 = smooth(c1Contour{ff,:},smoothPar);
    %sy2 = smooth(c2Contour{ff,:},smoothPar);
    sy2_1 = smooth(c2Contour{ff,:},2);
    sy1_1 = smooth(c1Contour{ff,:},2);

    [symax1,simax1,symin1,simin1] = extrema(sy1);
    [symax2,simax2,symin2,simin2] = extrema(sy2);
    [symax3,simax3,symin3,simin3] = extrema(sy2_1);
    [symax4,simax4,symin4,simin4] = extrema(sy1_1);

    sn = size(symax1,1)+size(symin1,1);
    sn2 = size(symax2,1)+size(symin2,1);
    sn2_1 = size(symax3,1)+size(symin3,1);
    M1_value(ff,2)=sn;
    M2_value(ff,2)=sn2;

    length = ymax1(1)-ymin1(1);
    width =  ymax2(1)-ymin2(1);
    perimeter = size(y1,1);

    M_features(ff,1)= length; %column1:length
    M_features(ff,2)= width;  %column2:width

    ROI_segment = maskofIMi(floor(ymin1:ymax1),floor(ymin2:ymax2));
    % figure, imshow(ROI_segment);
    nn = find(ROI_segment<1);
    area = size(nn,1);

    index = fix(symax4(1)); %??1
    contourX = contour(:,2);
    contourY = contour(:,1);

    No_simax = find(contour(:,1) == index); % meaning: simax
    plusfactor = contourX(max(No_simax))-contourX(min(No_simax));% ??2
    No_simax = find(contour(:,1) == index);

    startpoint = No_simax(max(find(No_simax < simax4(1))));
    endpoint = No_simax(min(find(No_simax > simax4(1))));

    while (plusfactor<=15);
        index = index -1;
        No_simax = find(contour(:,1) == index);

        startpoint = No_simax(max(find(No_simax < simax4(1))));
        endpoint = No_simax(min(find(No_simax > simax4(1))));

        plusfactor = contourX(endpoint)-contourX(startpoint);
    end

    if plusfactor >25
        index = index+1;

        No_simax = find(contour(:,1) == index);
        startpoint = No_simax(max(find(No_simax < simax4(1))));
        endpoint = No_simax(min(find(No_simax > simax4(1))));

        plusfactor = contourX(endpoint)-contourX(startpoint);
    end

%     endpoint = points(2);
%     startpoint = points(1);

    deltalength = abs(index - fix(symax1(1)));
    No_simax = find(contour == index);

%%%%%%%% replot the shape: start from the stem area %%%%%%%%%%%%%%%%%%%%%%%
    contourXX = cat(1,contourX(endpoint:size(contourX,1)),contourX(1:startpoint));
    contourYY = cat(1,contourY(endpoint:size(contourY,1)),contourY(1:startpoint));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% leaf tips
    contourX = contour(:,2);
    contourY = contour(:,1);
    %%%%%%%%%%%%%%%%%%%% edit 2015/06/04
    index_n = fix(simax3(1));
    no_idex_n = find(contour(:,1) ==  contour(index_n,1));
    
    if no_idex_n(1) < max(no_idex_n)
        contourStart = max(no_idex_n);
        contourend = min(no_idex_n);
        contour_of_leaf_x = cat(1,contourX(contourStart:end),contourX(1:contourend));
        contour_of_leaf_y = cat(1,contourY(contourStart:end),contourY(1:contourend));
        endi = size(contourY);
        H0 = 1:1:((endi(1)+contourend-contourStart)+1);
        H00 = 1:1:(endi(1)+contourend-contourStart);
        H1 = 1:((endi(1)+contourend-contourStart)+1)/400:abs((endi(1)+contourend-contourStart)+1);
        
        impoly_contour_of_leaf_x = interp1(H0,contour_of_leaf_x,H1);
        impoly_contour_of_leaf_y = interp1(H0,contour_of_leaf_y,H1);
        
    else
        contour_of_leaf_x = cat(1,contourX(contourStart:contourend));
        contour_of_leaf_y = cat(1,contourY(contourStart:contourend));
        
        H0 = 1:1:((contourend-contourStart)+1);
        H00 = 1:1:abs(contourend-contourStart);
        H1 = 1:(abs(contourend-contourStart)+1)/400:abs(abs(contourend-contourStart)+1);

        impoly_contour_of_leaf_x = interp1(H0,contour_of_leaf_x,H1);
        impoly_contour_of_leaf_y = interp1(H0,contour_of_leaf_y,H1);
    end
    
    clear Smooth_icolx
    clear Smooth_icolx
    Smooth_icolx(:,1) = smooth([impoly_contour_of_leaf_x]);
    Smooth_icolx(:,2) = smooth([impoly_contour_of_leaf_y]);
    slope = diff(Smooth_icolx);
    der_slope = diff(slope(:,end));
% 
% find the spectrum of ouline (tip) ---------------------------------------
    clear ALL
    ALL(:,2) = smooth([impoly_contour_of_leaf_y]);
    slopeALL = diff(ALL);
    der_slope_ALL = diff(slopeALL(:,end));

    Matrix_der_slope_all{ff} = der_slope_ALL;

% find the whole lead spectrum : from the bottom up : ???????????
    H0 = 1:1:(size(contourXX,1));
    H1 = 1:(size(contourXX,1))/160:(size(contourXX,1));
    clear impoly_contourX
    clear impoly_contourY
    clear ALLcontour
    impoly_contourX = interp1(H0,contourXX,H1);
    impoly_contourY = interp1(H0,contourYY,H1);
    ALLcontour(:,1) = smooth([impoly_contourX],10);
    ALLcontour(:,2) = smooth([impoly_contourY]);
    slopeofcontour = diff(ALLcontour);
    der_slopeofcontour1 = smooth(diff(slopeofcontour(:,1)),5);
    der_slopeofcontour2 = smooth(diff(slopeofcontour(:,2)),5);
    combine_der_slopeofcontour = der_slopeofcontour1+der_slopeofcontour2;

    % % 

% save matrix
Cal_features(ff,1) = sn;
Cal_features(ff,2) = sn2_1;
Cal_features(ff,3) = sn2;
Cal_features(ff,4) = perimeter;
Cal_features(ff,5) = area;
Cal_features(ff,6) = perimeter/area;
Cal_features(ff,7) = mean(mean(maskofIMi.*hsvIM4(:,:,3)));
% Cal_features(ff,8) = deltalength;
% Cal_features(ff,9) = length-deltalength;
Cal_features(ff,8) = (length-deltalength)/width;

if size(der_slope,1) <260
    Cal_features(ff,9) = 0; % ??????????????
else
Cal_features(ff,9) = max(der_slope(160:260));
end

der_Cell{ff,1} = der_slopeofcontour1(:,1);
der_Cell{ff,2} = der_slopeofcontour2(:,1);
der_Cell{ff,3} = combine_der_slopeofcontour(:,1);

%clear AAA
%toc

