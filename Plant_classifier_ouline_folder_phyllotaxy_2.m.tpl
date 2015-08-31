%% ^^^^^ configuration zone

%% Run the file 
tic % count (timer)
%% prev edition script legacy var ff
ff = 1; % read the specific file

cd(exec_dir)
dir_list=dir([upload_dir '/**'] );
%dir_list=dir([exec_dir '/' dir_list(dd).name '/IM***'])

    % select the subfolder and run the data 
   % lock the line 9-17
   cd(exec_dir)
   run Leaf_Feature_extract
        
   % cd(exec_dir) 
   load('1_617.mat') % training data
    class = textread('ATextName.txt','%s','delimiter','\t'); % Specify each training data (species)
    
    
    %% Calculate the ouline

    Sum_sqerror_1 = zeros(size(Cal_features_sample,1),1); % declare a space for calculating least square error : dimension 1  
    Sum_sqerror_2 = zeros(size(Cal_features_sample,1),1); % declare a space for calculating least square error : dimension 2
    Sum_sqerror_3 = zeros(size(Cal_features_sample,1),1); % declare a space for calculating least square error : dimension 3
    
    % Sum_squerror_1 = size(trainingData, filesSubfolder)
    for i = 1:size(Cal_features_sample,1) 
        for ii = 1:size(der_Cell,1) 
            dut_1 = der_Cell{ii,1}; % extract the LSQR data from training data: dimension 1
            dut_2 = der_Cell{ii,2}; % extract the LSQR data from training data: dimension 2
            dut_3 = der_Cell{ii,3}; % extract the LSQR data from training data: dimension 3

            sample_1 = der_Cell_sample{i,1}; % test data: dimension 1 
            sample_2 = der_Cell_sample{i,2}; % test data: dimension 2
            sample_3 = der_Cell_sample{i,3}; % test data: dimension 3

            Sum_sqerror_1(i,ii) = sum((dut_1(1:150)-sample_1(1:150)).^2); % calculate LSQR: dimension 1
            Sum_sqerror_2(i,ii) = sum((dut_2(1:150)-sample_2(1:150)).^2); % calculate LSQR: dimension 2
            Sum_sqerror_3(i,ii) = sum((dut_3(1:150)-sample_3(1:150)).^2); % calcualte LSQR: dimension 3
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    %% species: count how many spicies are there in the training data
    % ----------------------------------------------------------------------------------------------------------------------
        clear a b k ansTest 
        k = 1; % count = 1
        b = zeros(size(Cal_features_sample,1),1);
        a = zeros(size(Cal_features_sample,1),1);
    
        for i = 1:size(Cal_features_sample,1)   
            if b(i) == 0;   % accumulate each judge
                a(strcmp(class, class{i}))=k; % per time judge
                b(find(a==k))=k; 
                k = k+1;
                
                %Cell_specis{k-1,1} = class{i};
                %Cell_phyllotaxy{k,2} = Cal_features_sample(i,1)
                ALLtypephyllotaxy(k-1,1) = Cal_features_sample(i,1);
            end
        end
    
    % ----------------------------------------------------------------------------------------------------------------------
    %% phyllotaxy: based on the information of phyllotaxy, delete the unrelated species
    % ----------------------------------------------------------------------------------------------------------------------
    
    if ceil(mean(Cal_features(:,1))) > 5
       retainData = [1, 2, 15, 19, 34, 57, 59, 61, 62]'; 
       deleteData = [1:max(b)]'; deleteData(retainData,:) =[];
    else 
       deleteData = find(ALLtypephyllotaxy ~= phyllotaxy);
       retainData = find(ALLtypephyllotaxy == phyllotaxy);
    end
    
    %%
    SpeciesK_1 = zeros(size(Sum_sqerror_1,2),max(b));
    SpeciesK_2 = zeros(size(Sum_sqerror_1,2),max(b));
    SpeciesK_3 = zeros(size(Sum_sqerror_1,2),max(b));

    for ii = 1:size(Sum_sqerror_1,2)
        for i = 1:max(b)
            trail_1 = Sum_sqerror_1(:,ii);
            trail_2 = Sum_sqerror_2(:,ii);
            trail_3 = Sum_sqerror_3(:,ii);

            SpeciesK_1(ii,i) = mean(trail_1(find(b==i)));
            SpeciesK_2(ii,i) = mean(trail_2(find(b==i)));
            SpeciesK_3(ii,i) = mean(trail_3(find(b==i)));   
        %     resultresult{ii,1} = class{min(find(b==r(1)))}
        end
    end

    SpeciesK_1(:,deleteData) = [];
    SpeciesK_2(:,deleteData) = []; 
    SpeciesK_3(:,deleteData) = [];
    
    clear r_matrix_1 r_matrix_2 r_matrix_3
    clear c_matrix_1 c_matrix_2 c_matrix_3
    
    r = zeros(size(Sum_sqerror_1,2),max(b));
    for ii = 1:size(Sum_sqerror_1,2)
        [c1 r1] = sort(SpeciesK_1(ii,:),'ascend'); %c1 = value, r1 = serial number
        [c2 r2] = sort(SpeciesK_2(ii,:),'ascend'); %c1 = value, r1 = serial number
        [c3 r3] = sort(SpeciesK_3(ii,:),'ascend'); %c1 = value, r1 = serial number  
        
        r_matrix_1(ii,:)= retainData(r1(:));
        r_matrix_2(ii,:)= retainData(r2(:));
        r_matrix_3(ii,:)= retainData(r3(:));
        
        c_matrix_1(ii,:)= c1;
        c_matrix_2(ii,:)= c2;
        c_matrix_3(ii,:)= c3;

    end
    
   
    clear resultresult1 resultresult2 resultresult3
    
   
    if size(r_matrix_1,2) < 10
        rank = size(r_matrix_1,2);
    else rank = 10;
    end
    
    for ii=1:rank %rank ?__?
        for jj = 1:size(Sum_sqerror_1,2)    %number of files
        resultresult1{jj,ii} = class{min(find(b==r_matrix_1(jj,ii)))};
        resultresult2{jj,ii} = class{min(find(b==r_matrix_2(jj,ii)))};
        resultresult3{jj,ii} = class{min(find(b==r_matrix_3(jj,ii)))};
        
        resultA_1(jj,ii) = r_matrix_1(jj,ii);   % species ?? #
        resultA_2(jj,ii) = r_matrix_2(jj,ii);
        resultA_3(jj,ii) = r_matrix_3(jj,ii);        

        end
    end
    
    for ii=1:rank %rank ?__?
        for jj = 1:size(Sum_sqerror_1,2)
        resultrmatrix1{jj,ii} = c_matrix_1(jj,ii);
        resultrmatrix2{jj,ii} = c_matrix_2(jj,ii);
        resultrmatrix3{jj,ii} = c_matrix_3(jj,ii);
        end
    end

    
    % ----------------------------------------------------------------------------------------------------------------------
    %% Bayes classification: 
    % ----------------------------------------------------------------------------------------------------------------------
    load('1_617.mat');
    k = 1;
    class = textread('ATextName.txt','%s','delimiter','\t');

    attrib = Cal_features_sample';         % Standard Sample
    pcata = zeros(max(b(:)),1);     % size = type of plants
    pcata(1:end,1)=1/max(b(:));

    p=zeros(max(b(:)),1);
    
    insertC = ones(size(Cal_features,1),1)*phyllotaxy;
    Cal_features_2 = cat(2,insertC,Cal_features);
    %the data we need test:
    %%
    clear test_data
    for ii = 1:size(Cal_features_2,1)
        test_data = Cal_features_2(ii,:);
        for i = 1:size(p,1)
            allofonecata=0;
           cp=1;
           col=1;
           coli=1;
           meet=zeros(1,size(attrib,2));
            for j=1:size(attrib,2)
               %compute the Conditional probability
               if a(j,1)==i
                      allofonecata=allofonecata+1;
                      for coli=1:size(attrib,1)
                         if attrib(coli,j)==test_data(1,coli)
                             meet(1,coli)=meet(1,coli)+1;
                         end
                      end
               end
            end
            for coli=1:size(attrib,1)
             cp=cp*((meet(1,coli)/allofonecata)+0.01);
            end
            p(i,1)=cp*pcata(i,1);
        end
        Result(:,ii) = p(:);
    end

    Result(deleteData,:) = [];

    %?????????
    max=p(1,1);
    resulti=1;
    for ii = 1:size(Cal_features,1);
        for i=1:size(p,1);
            [c r] = sort(Result(:,ii),'descend');
            resultrank = retainData(r(:));
        end    
            if size(Result,1) < 10;
                rank = size(r_matrix_1,2);
                else rank = 10;
            end
            for j = 1:rank;
            resultresult{ii,j} = class{min(find(b==resultrank(j)))}; % ??
            result_bayes_cvalue{ii,j} = c(j)*10^15;         % ??
            resultB_1(ii,j) = resultrank(j);
            end
        
    end

    % ----------------------------------------------------------------------------------------------------------------- 
    %% CALCULATE ALL
    % -----------------------------------------------------------------------------------------------------------------
    clear acumulatorSpecies
    for ii = 1:size(Cal_features,1);
        for i=1:size(p,1) ;   
               acumulatorSpecies(i) = mean([find(resultA_1(ii,:)==i) ...
                  find(resultA_2(ii,:)==i)...
                  find(resultA_3(ii,:)==i)...
                  find(resultB_1(ii,:)==i)]);
        end
       [c r] = sort(acumulatorSpecies,'ascend');
       species(ii,:)= r(1:rank);
    end
    species

    %% create text file tell who is it
    cd(upload_dir)

	    header = 'Sample';

	    
	    MAT = [num2str(species(1,:))];
	    
	    ss = num2str(ff);
	    nur=[repmat(['0'],1,4-size(ss,2)),ss];
	    str = [header,'_',nur];
	    name = [str '.txt'];

	    cd(upload_dir);
	    dlmwrite(name, MAT(1:end),'delimiter','%s\t','newline','pc');

toc % timer end
