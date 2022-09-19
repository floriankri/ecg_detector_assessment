function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
    % By Jan Gieraltowski, Piotr Podziemski
    % ---- check size of ECG - do the input from csv is horizontal or vertical? ----
    if size(ECG,2)>size(ECG,1)
        ECG = ECG';
    end
    %------------PARAMETRY--------------------------------------------
    fs      = 1000;              % sampling frequency
    N       = size(ECG,2);       % number of abdominal channels
    Nlength       = size(ECG,1); % length of the recording (in samples)
    
    %------------PARAMETRY DO OPTYMALIZAZJI---------------------------    
    window=100;                  % length of the moving median window (in samples); Moving median is used for the de-trending purposes.
    zakresIntersecta=2;          % used when choosing two good channels based on an intersection. Two annotations treated as the same cannot differ more than 3 ms.
    lowLength=7;                 % minimal length of the fECG R-peak repolarization
    upLength=20;                 % maximal length of the fECG R-peak repolarization
    czyZero = 0;                 % boolean; 1: check, if repolarization is below zero - level during search for peaks
    rangeSD = 700;               % maximum limit of SD that can occur in the skeleton set of RR peaks
    rangeSD_20 = 1000;
    skeletonLimit= 75;           % maximum deviation of the interval from the mean interval in the chosen skeleton (only during first BRUTFORS stage)
    oknoDoboruPiku= 100;         % maximum deviation from the mean value of the interval at the border of a single skeleton piece
    lowZakresMODA=300;           % minimal value of the good mean value of intervals in the skeleton(s) (used at almost all stages)
    upZakresMODA=525;            % maximum value of the good mean value of intervals in the skeleton(s) (used at almost all stages)
    dziuraLimit = 600;           % minimal length of a hole in skeleton to be fit during last stage    
    %--1.WYPELNIJ PUSTKI---------------------------------------------- 
    ECGstage0= wypelniaczPustek(tm,ECG);            % fiting the empty potential trace with polynomial of order 2
    %ECGstageTEMP = movingMedian(ECGstage0,window);
    %--2.Fourier Notch filter ----------------------------------------
    NFFT = 2^nextpow2(Nlength);
    f = 1000/2*linspace(0,1,NFFT/2+1);
    f=f';
    for iterchannel =1:4
        Y = fft(ECGstage0(:,iterchannel),NFFT)/Nlength; 
        FURIER = [f, abs(Y(1:NFFT/2+1)).*abs(Y(1:NFFT/2+1))];        
        %plot(FURIER(:,1),FURIER(:,2));
        %pause;
        [~,pos49] = min(abs(f-49));
        [~,pos51] = min(abs(f-51));
        [~,pos59] = min(abs(f-59));
        [~,pos61] = min(abs(f-61));
        [~,pos47] = min(abs(f-47));
        [max50Hz,pos_max50Hz] = max(FURIER(pos49:pos51,2));
        pos_max50Hz=pos_max50Hz+pos49;
        [max60Hz,pos_max60Hz] = max(FURIER(pos59:pos61,2));
        pos_max60Hz=pos_max60Hz+pos59;
         
        porownanie = max(FURIER(pos47:pos49,2));
        if( max50Hz > porownanie*30)
            ECGstage0(:,iterchannel)= notchFilter(ECGstage0(:,iterchannel),50);  
           
        elseif( max60Hz > porownanie*30)
            ECGstage0(:,iterchannel)= notchFilter(ECGstage0(:,iterchannel),60);
           
        end
    end  
    
    %--2.RUCHOMAMEDIANA-----------------------------------------------
    ECGstage1 = movingMedian(ECGstage0,window);

    
    %--3.PIERWSZY BRUTFORS-------------------------------------------- 
    celsOfBest = cell(1000,8);                      % cell array of the results of each brutfors loop
    counterek=1;                                    % counter of all the loops of parameter space chceck
    %brutMatrix = ECGstage1;                        
    brutMatrix = [ECGstage1, -ECGstage1];           % matrix of all signals (normal and revese polarisation)
    %size(brutMatrix)           % matrix of all signals (normal and revese polarisation)
    pikiBrutfors1 = cell(1,size(brutMatrix,2));     % cell array for detected RR peaks for 8 signals (size(brutMatrix,2) = 8!)
    permtable = nchoosek(1:size(brutMatrix,2),2);   % table of possible permutation of 8 signals.
    tabelaDlugosci =zeros(size(permtable,1),2);     % table of the number of annotation for one pair from all permutations,
    %tabelaBrutforsow = cell(size(permtable,1),2); 
    
    %--4.Main loop----------------------------------------------------
    % dolnagranica - the lower percentile of the distribution of recorded
    %                potential values, in which we search for deceleration
    %                start.
    % gornagranica - the upper percentile of the distribution of recorded
    %                potential values, in which we search for deceleration
    %                start. Cannot be lower than 0.9
    %
    for dolnagranica=0.99:-0.02:0.79
        startowa=dolnagranica+0.01;
        if startowa < 0.90
            startowa = 0.90;
        end
        for gornagranica = 0.99:-0.02:startowa
            % here we search for peaks in all of the 8 channels. We store
            % the result in pikiBrutfors1
            for itertemp=1:size(brutMatrix,2)
                pikitemp=findDeccelerations(brutMatrix(:,itertemp),dolnagranica,gornagranica,lowLength,upLength,czyZero);
                pikiBrutfors1(1,itertemp) = {pikitemp};
            end
            % now we search for the biggest intersection between two signals
            for itertemp=1:size(permtable,1)
                A = cell2mat(pikiBrutfors1(permtable(itertemp,1)));
                B = cell2mat(pikiBrutfors1(permtable(itertemp,2)));
                C = intersect(A,B);
                for iterP=1:zakresIntersecta
                    if(isempty(C))
                        C=[];
                    end
                    temp1=intersect(A,B+iterP);
                    if(isempty(temp1))
                    else
                     C = cat(1,C,temp1);
                    end
                    temp2=intersect(A,B-iterP);
                    if(isempty(temp2))
                    else
                         C = cat(1,C,temp2);
                    end
                end
                tabelaDlugosci(itertemp,1)=itertemp;
                tabelaDlugosci(itertemp,2)=size(C,1); 
            end
            [~,iii]=max(tabelaDlugosci(:,2));
                   
            % here, two mixed signals from chosen permutation will be temporary stored
            ECGstageBrut1 = (brutMatrix(:,permtable(iii,1)) + brutMatrix(:,permtable(iii,2)))/2;
            % here, first signal from chosen permutation will be temporary stored    
            ECGstageBrut2 = (brutMatrix(:,permtable(iii,1)));
            % here, second signal from chosen permutation will be temporary stored
            ECGstageBrut3 = (brutMatrix(:,permtable(iii,2)));

            pikifound1 = findDeccelerations(ECGstageBrut1,dolnagranica,gornagranica,lowLength,upLength,czyZero);
            pikifound2 = findDeccelerations(ECGstageBrut2,dolnagranica,gornagranica,lowLength,upLength,czyZero);
            pikifound3 = findDeccelerations(ECGstageBrut3,dolnagranica,gornagranica,lowLength,upLength,czyZero);

            celsOfBest(counterek,1)={dolnagranica};
            celsOfBest(counterek,2)={gornagranica};
            celsOfBest(counterek,3)={pikifound1};
            celsOfBest(counterek,4)={pikifound2};
            celsOfBest(counterek,5)={pikifound3};
            celsOfBest(counterek,6)={permtable(iii,1)};
            celsOfBest(counterek,7)={permtable(iii,2)};
            celsOfBest(counterek,8)={iii};
            counterek=counterek+1;
        end
    end
    %--5.BRUTFORS-POROWNANIE    
    tablicaSzkieletow=zeros(size(celsOfBest(:,1),1),3);
    for iterGlowny=1:size(celsOfBest(:,1))
        for iterKanal=1:3
            A = cell2mat(celsOfBest(iterGlowny,2+iterKanal));
            MODA = diff(A);
            if(size(MODA,1)>=1)
                MODA(MODA>upZakresMODA) =[];
                MODA(MODA<lowZakresMODA) =[];
            end
            MEANMODA = mean(MODA);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ANNour = A;
            dANNour =diff(A);
            if(size(dANNour,1)>1)
                for iterS=2:size(dANNour,1)
                    if dANNour(iterS) < MEANMODA+skeletonLimit && dANNour(iterS) > MEANMODA-skeletonLimit &&...
                       dANNour(iterS-1) < MEANMODA+skeletonLimit && dANNour(iterS-1) > MEANMODA-skeletonLimit
                    else
                        ANNour(iterS)=999999;
                    end
                end
            end
            ANNour(ANNour==999999)=[];
            ANNour = sort(ANNour);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(isempty(ANNour))
                tablicaSzkieletow(iterGlowny,iterKanal)=0;
            else
                tablicaSzkieletow(iterGlowny,iterKanal)=size(ANNour,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %----warnuek na SD-------
            if std(dANNour) > rangeSD
                tablicaSzkieletow(iterGlowny,iterKanal) = NaN;
            end      
        end
    end
    [CCC,I] = max(tablicaSzkieletow);
    [~,whichchannel] = max(CCC);
    %--5.BRUTFORS-UZUPELNIENIA SZKIELETU   
    tablicaSzkieletowDoUzupelnien=zeros(size(tablicaSzkieletow,1),4);
    licznik=1;
    for iterSzkielety=1:size(tablicaSzkieletow,1)    
        if (tablicaSzkieletow(iterSzkielety,1)>0.3*max(CCC) )           
                tablicaSzkieletowDoUzupelnien(licznik,1)=cell2mat(celsOfBest(iterSzkielety,8)); % ktora permutacja
                tablicaSzkieletowDoUzupelnien(licznik,2)=cell2mat(celsOfBest(iterSzkielety,1)); % dolna granica
                tablicaSzkieletowDoUzupelnien(licznik,3)=cell2mat(celsOfBest(iterSzkielety,2)); % gorna granica
                tablicaSzkieletowDoUzupelnien(licznik,4)=1; % suma
                licznik=licznik+1;
        end
        if (tablicaSzkieletow(iterSzkielety,2)>0.3*max(CCC))           
                tablicaSzkieletowDoUzupelnien(licznik,1)=cell2mat(celsOfBest(iterSzkielety,8)); % ktora permutacja
                tablicaSzkieletowDoUzupelnien(licznik,2)=cell2mat(celsOfBest(iterSzkielety,1)); % dolna granica
                tablicaSzkieletowDoUzupelnien(licznik,3)=cell2mat(celsOfBest(iterSzkielety,2)); % gorna granica
                tablicaSzkieletowDoUzupelnien(licznik,4)=2; % 1 kanal
                licznik=licznik+1;
        end
        if (tablicaSzkieletow(iterSzkielety,3)>0.3*max(CCC))        
                tablicaSzkieletowDoUzupelnien(licznik,1)=cell2mat(celsOfBest(iterSzkielety,8)); % ktora permutacja
                tablicaSzkieletowDoUzupelnien(licznik,2)=cell2mat(celsOfBest(iterSzkielety,1)); % dolna granica
                tablicaSzkieletowDoUzupelnien(licznik,3)=cell2mat(celsOfBest(iterSzkielety,2)); % gorna granica
                tablicaSzkieletowDoUzupelnien(licznik,4)=3; % 2 kanal
                licznik=licznik+1;
        end
    end
    
    lowg = cell2mat(celsOfBest(I(1,whichchannel),1));
    highg = cell2mat(celsOfBest(I(1,whichchannel),2));
    ktorapermutacja = cell2mat(celsOfBest(I(1,whichchannel),8));
    %--2.CHANNELMIXING-------------------------    
    ECGstage2 = (brutMatrix(:,permtable(ktorapermutacja,1)) + brutMatrix(:,permtable(ktorapermutacja,2)))/2;
    sig1 = (brutMatrix(:,permtable(ktorapermutacja,1)));
    sig2 = (brutMatrix(:,permtable(ktorapermutacja,2)));
    
    if whichchannel==1
        dobrySYGNAL=ECGstage2;
    elseif whichchannel==2
        dobrySYGNAL=sig1;
    else
        dobrySYGNAL=sig2;
    end
     
    %%%BRUTFORS2-------------------------------------------------------------------
    celsOfBest2 = cell(1000,4);
    counterek=1; 
    for dolnadlugosc=5:1:12
        for gornadlugosc = 15:1:30
            ktorezero =0.0;
            pikifound = findDeccelerations(dobrySYGNAL,lowg,highg,dolnadlugosc,gornadlugosc,ktorezero);
            celsOfBest2(counterek,1)={dolnadlugosc};
            celsOfBest2(counterek,2)={gornadlugosc};
            celsOfBest2(counterek,3)={ktorezero};
            celsOfBest2(counterek,4)={pikifound};
            counterek=counterek+1;

            ktorezero = 1.0;
            pikifound = findDeccelerations(dobrySYGNAL,lowg,highg,dolnadlugosc,gornadlugosc,ktorezero);
            celsOfBest2(counterek,1)={dolnadlugosc};
            celsOfBest2(counterek,2)={gornadlugosc};
            celsOfBest2(counterek,3)={ktorezero};
            celsOfBest2(counterek,4)={pikifound};
            counterek=counterek+1;
        end
    end
    %--3.BRUTFORS-POROWNANIE dlugosci
    tabliceSzkieletow2=zeros(size(celsOfBest2(:,1),1),1);
    for iter=1:size(celsOfBest2(:,1))
        A = cell2mat(celsOfBest2(iter,4));
        MODA = diff(A);
        if(size(MODA,1)>=1)
            MODA(MODA>upZakresMODA) =[];
            MODA(MODA<lowZakresMODA) =[];
        end        
        if(isempty(MODA))
            tabliceSzkieletow2(iter)=0;
        else
            tabliceSzkieletow2(iter)=size(MODA,1);
        end
        MEANMODA = mean(MODA);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ANNour = A;
        dANNour =diff(A);
        if(size(dANNour,1)>1)
            for iterS=2:size(dANNour,1)
                %XXXX zacina sie przy dlugosci = 1
                if dANNour(iterS) < MEANMODA+skeletonLimit && dANNour(iterS) > MEANMODA-skeletonLimit &&...
                   dANNour(iterS-1) < MEANMODA+skeletonLimit && dANNour(iterS-1) > MEANMODA-skeletonLimit
                else
                    ANNour(iterS)=999999;
                end
            end
        end
        ANNour(ANNour==999999)=[];
        ANNour = sort(ANNour);
        if(isempty(ANNour))
            tabliceSzkieletow2(iter)=0;
        else
            tabliceSzkieletow2(iter)=size(ANNour,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----warnuek na SD-----------
        if std(dANNour) > rangeSD
            tabliceSzkieletow2(iter) = NaN;
        end         
    end    
    [~,I] = max(tabliceSzkieletow2);    
    lowLengthBest = cell2mat(celsOfBest2(I(1,1),1));
    highLengthBest = cell2mat(celsOfBest2(I(1,1),2));
    zeroBest = cell2mat(celsOfBest2(I(1,1),3));
    
    ANNourREF = findDeccelerations(dobrySYGNAL,lowg,highg,lowLengthBest,highLengthBest,zeroBest);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%MODYFIKACJA SIERPNIOWA 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [motherAnn, motherECG] = usuwaczMatki( dobrySYGNAL);
    odjetaMatka = dobrySYGNAL-motherECG;   
       
    przefiltrowany= notchFilter(dobrySYGNAL,50); 
    signalZeroMother = zeroMother(motherAnn, przefiltrowany, -140, 280);
    indexMax=1;
    [~,indexMax]=max(abs(signalZeroMother(1000:size(signalZeroMother,1)-1000)));
    

    
    
    liczKorelacje = 1;
	if size(ANNourREF,1) < 3
		ANNourREF = findDeccelerations(signalZeroMother,0.9,1.0,5,30,zeroBest);
		ANNourREF2 = findDeccelerations(-signalZeroMother,0.9,1.0,5,30,zeroBest);
        cleantemp1=annCleaning(ANNourREF, 7000,lowZakresMODA,7000,7000);
        cleantemp2=annCleaning(ANNourREF2, 7000,lowZakresMODA,7000,7000);
        if(size(cleantemp1,1) < size(cleantemp2,1))
            cleantemp1 = cleantemp2;
            ANNourREF = ANNourREF2;
            signalZeroMother=-signalZeroMother;
        end
        ANNour2=cleantemp1;
        liczKorelacje = 0;
    else
        ANNour2=annCleaning(ANNourREF, upZakresMODA,lowZakresMODA,oknoDoboruPiku,skeletonLimit);       
    end
    
 	if size(ANNour2,1) < 3
        indexMax=indexMax+1000;
        ANNourRESCUE=[];
        for iterM=indexMax:425:60000
            ANNourRESCUE(end+1,1) = iterM;
        end
        for iterM=indexMax:-425:1
            ANNourRESCUE(end+1,1) = iterM;
        end
        ANNourRESCUE=sort(ANNourRESCUE);
        ANNour2=ANNourRESCUE;
	end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % KONIEC - MAMY NAJLEPIEJ ZNALEZIONE PIKI
    % teraz - dobor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    MEANMODA = 425;   
    ANNourSTAGE2 = ANNour2;
    
    % korelacja template'u dziecka plus template do QT
    korelacjaDziecko=templateOfChild(odjetaMatka, ANNourSTAGE2);
    poSkorelowaniuZTemplatemDziecka = abs(korelacjaDziecko).*odjetaMatka;
    

    celsOfBest3 = cell(1000,3);
    counterek=1; 
    for dolnagranica=0.99:-0.01:0.80
        startowa=dolnagranica+0.01;
        if startowa < 0.90
            startowa = 0.90;
        end
        for gornagranica = startowa:0.01:1.0
                    pikifound = findDeccelerations(poSkorelowaniuZTemplatemDziecka,dolnagranica,gornagranica,7,20,10000);
                    celsOfBest3(counterek,1)={dolnagranica};
                    celsOfBest3(counterek,2)={gornagranica};
                    celsOfBest3(counterek,3)={pikifound};
                    counterek=counterek+1;
        end
    end
    %--3.BRUTFORS-POROWNANIE 
    [firstBest, secondBest] = brutforsPorownanie(celsOfBest3, upZakresMODA,lowZakresMODA,rangeSD_20,skeletonLimit);  
    lowBestPoUsunieciuMatki = firstBest;  
    upBestPoUsunieciuMatki = secondBest;
    
    %2 brutfors dlugosci
    celsOfBest4 = cell(1000,3);
    counterek=1; 
    for dolnadlugosc=5:1:12
        for gornadlugosc = 15:1:30
                    pikifound = findDeccelerations(poSkorelowaniuZTemplatemDziecka,lowBestPoUsunieciuMatki,upBestPoUsunieciuMatki,dolnadlugosc,gornadlugosc,10000);
                    celsOfBest4(counterek,1)={dolnadlugosc};
                    celsOfBest4(counterek,2)={gornadlugosc};
                    celsOfBest4(counterek,3)={pikifound};
                    counterek=counterek+1;
        end
    end
    %--3.BRUTFORS-POROWNANIE 
    [firstBest, secondBest] = brutforsPorownanie(celsOfBest4, upZakresMODA,lowZakresMODA,rangeSD_20,skeletonLimit);  
    lowLengthBestPoUsunieciuMatki = firstBest;  
    upLengthBestPoUsunieciuMatki = secondBest;
    
    % szkielet 2.0
    ANNourREF20 = findDeccelerations(poSkorelowaniuZTemplatemDziecka,lowBestPoUsunieciuMatki,upBestPoUsunieciuMatki,lowLengthBestPoUsunieciuMatki,upLengthBestPoUsunieciuMatki,10000);
    
	if(isempty(ANNourREF20))
		ANNourREF20 = findDeccelerations(poSkorelowaniuZTemplatemDziecka,0.6,0.99,5,40,0);
    end   
	if(size(ANNourREF20,1) < 3)
		ANNourREF20 = findDeccelerations(poSkorelowaniuZTemplatemDziecka,0.6,0.99,5,40,0);
    end
    ANNour_s20granice=annCleaning(ANNourREF20, upZakresMODA,lowZakresMODA,oknoDoboruPiku,skeletonLimit);
    ANNour_szkielt20 = ANNour_s20granice;
    dANNour_szkielt20=diff(ANNour_s20granice);

    % SUMA SZKIELETOW
    if(liczKorelacje~=0)
        polaczonySzkielet=skeletonJoin(ANNour2,ANNour_szkielt20,lowZakresMODA);
    else
        polaczonySzkielet=ANNour2;
    end
    
    %polaczonySzkielet=annCleaning(polaczonySzkielet, upZakresMODA,lowZakresMODA,oknoDoboruPiku,skeletonLimit);

    
    %disp('odejmowanie matek');
    %tic
    brutMatrixBezMatki = brutMatrix;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%SUMOWANIE SZKIELETOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tablicaSzkieletowDoUzupelnien(tablicaSzkieletowDoUzupelnien==0)=[];
    tablicaSzkieletowDoUzupelnien = reshape(tablicaSzkieletowDoUzupelnien,[],4);    
    szkieletyDoDodania = cell(size(tablicaSzkieletowDoUzupelnien,1),1);                      % cell array of the results of each brutfors loop
    kolenoscSzkieletow = zeros(size(tablicaSzkieletowDoUzupelnien,1),2);
    %tic
    if(size(tablicaSzkieletowDoUzupelnien,1)>1)
        for licznik=1:size(tablicaSzkieletowDoUzupelnien,1)

             lowgtemp = tablicaSzkieletowDoUzupelnien( licznik,2);
             highgtemp = tablicaSzkieletowDoUzupelnien( licznik,3);
             ktorapermutacjatemp = tablicaSzkieletowDoUzupelnien( licznik,1);
             kanaltemp  = tablicaSzkieletowDoUzupelnien( licznik,4);

            if kanaltemp==1
                testSYGNAL=(brutMatrixBezMatki(:,permtable(ktorapermutacjatemp,1)) + brutMatrixBezMatki(:,permtable(ktorapermutacjatemp,2)))/2;
            elseif kanaltemp==2
                testSYGNAL=(brutMatrixBezMatki(:,permtable(ktorapermutacjatemp,1)));
            else
                testSYGNAL=(brutMatrixBezMatki(:,permtable(ktorapermutacjatemp,2)));
            end
            
             pikitemp= findDeccelerations( testSYGNAL,lowgtemp,highgtemp,lowLengthBest,highLengthBest,zeroBest);
             pikiclean=annCleaning(pikitemp, upZakresMODA,lowZakresMODA,oknoDoboruPiku,skeletonLimit);

             szkieletyDoDodania(licznik,1)={pikiclean};
             kolenoscSzkieletow(licznik,1)=licznik;
             kolenoscSzkieletow(licznik,2)=size(pikiclean,1);
        end
    end
    kolenoscSzkieletow=sortrows(kolenoscSzkieletow,-2);
    tempSzkielet=polaczonySzkielet;
    %kolenoscSzkieletow
    if(size(kolenoscSzkieletow,1)>1)
        for licznik=1:size(kolenoscSzkieletow,1)
                tempSzkielet=skeletonJoin(tempSzkielet,cell2mat(szkieletyDoDodania(kolenoscSzkieletow(licznik,1))),lowZakresMODA);
        end
    end
    polaczonySzkielet=tempSzkielet;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [templateQT,ourQT,positionOfQ,positionOfT] =estimateQT(odjetaMatka, polaczonySzkielet);
    ANNourSTAGE2 = polaczonySzkielet;
    dANNourSTAGE2=diff(polaczonySzkielet);
    
    
    
    
%       etap czwarty - usun odstajace
    MEANMODAdoSTD = 425;
    MEANSTD = 30;
    MODA = diff(ANNourSTAGE2);
    if(size(MODA,1)>=1)
        MODA(MODA>upZakresMODA) =[];
        MODA(MODA<lowZakresMODA) =[];
        MEANMODAdoSTD=mean(MODA);
        MEANSTD = 1.5*std(MODA);
    end       
    
    for iterSH=2:size(dANNourSTAGE2)-2
           
         MODAITER = diff(ANNourSTAGE2);
         if (size(MODAITER,1)>=1)
            MODAITER_TABLE=[];
            licznik_d=1;
            iter_d=iterSH-1; 
            if iter_d < 1
                iter_d=1;
            end
            while licznik_d<=5;
                if(MODAITER(iter_d) <upZakresMODA && MODAITER(iter_d) >lowZakresMODA )
                    MODAITER_TABLE(end+1,1) = MODAITER(iter_d);
                    licznik_d=licznik_d+1;
                end                                    
                iter_d=iter_d-1;  
                if iter_d < 1
                    break;
                end
            end
            licznik_g=1;
            iter_g=iterSH+1; 
            if iter_g >= size(MODAITER,1)
            	iter_g=size(MODAITER,1);
            end
            while licznik_g<=5;
                if(MODAITER(iter_g) <upZakresMODA && MODAITER(iter_g) >lowZakresMODA )
                    MODAITER_TABLE(end+1,1) = MODAITER(iter_g);
                    licznik_g=licznik_g+1;
                end   
                iter_g=iter_g+1;  
                if iter_g >= size(MODAITER,1)-1
                     break;
                end
            end                                
            if(isempty(MODAITER_TABLE))
                MODAITER_TABLE = MEANMODAdoSTD;
            end  
            %XXXXXXXXXXXXXXXXXXx
            MEANMODAdoSTD=mean(MODAITER_TABLE);
        end   
    
        if dANNourSTAGE2(iterSH) < dziuraLimit && ( dANNourSTAGE2(iterSH) >= MEANMODAdoSTD + MEANSTD || dANNourSTAGE2(iterSH) <= MEANMODAdoSTD - MEANSTD   )      
             if(dANNourSTAGE2(iterSH-1)>600)
                 ANNourSTAGE2(iterSH)=999999;
             else
                 ANNourSTAGE2(iterSH+1)=999999;
             end
        end
    end
    ANNourSTAGE2(ANNourSTAGE2==999999)=[];
    polaczonySzkielet = ANNourSTAGE2;
    dANNourSTAGE2=diff(ANNourSTAGE2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%MODYFIKACJA SIERPNIOWA KONIEC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %etap 4 - MAMY SZKIELET
    
    if(size(ANNourSTAGE2,1)>1 && size(ANNourSTAGE2,2)>0) 
        stopWatch = 1;
        while (stopWatch == 1)
            stopWatch =0; 
            fills=[];

            for iter=1:size(dANNourSTAGE2,1)
                 if iter == 1 && ANNourSTAGE2(iter)-1 > MEANMODA
                     fitTABLE=zeros(4,2);
                     if(iter-2>0)
                        fitTABLE(1,1)=ANNourSTAGE2(iter-2);
                        fitTABLE(1,2)=dANNourSTAGE2(iter-2);
                     end                         
                     if(iter-1>0)
                        fitTABLE(2,1)=ANNourSTAGE2(iter-1);
                        fitTABLE(2,2)=dANNourSTAGE2(iter-1);
                     end                         
                     if(iter+1<=size(dANNourSTAGE2,1))
                        fitTABLE(3,1)=ANNourSTAGE2(iter);
                         fitTABLE(3,2)=dANNourSTAGE2(iter);
                     end
                     if(iter+2<=size(dANNourSTAGE2,1))
                        fitTABLE(4,1)=ANNourSTAGE2(iter+1);
                        fitTABLE(4,2)=dANNourSTAGE2(iter+1);
                     end
                     fitTABLE(all(fitTABLE==0,2),:)=[];
                     warning('off','all');
                     p = polyfit(fitTABLE(:,1),fitTABLE(:,2),1);
                     warning('on','all'); 
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     aktualnaSZPILA=1;
                     nastepnaSZPILA=floor(ANNourSTAGE2(iter));
                     if nastepnaSZPILA>Nlength
                         nastepnaSZPILA=Nlength;
                     end    
                     wycinekSYGNALU=dobrySYGNAL(aktualnaSZPILA:nastepnaSZPILA,1);
                     luznePIKI=findDeccelerations(wycinekSYGNALU,0.8,1,5,40,10000000);                
                     luznePIKI=luznePIKI-1;
                     if(isempty(luznePIKI))    
                        fills = cat(1,fills,round(nastepnaSZPILA-mean(fitTABLE(:,2))));
                        stopWatch =1;
                     else
                         tablicaPRZECIEC=zeros(size(luznePIKI,1),1);
                         %%%%%%%%%%%%%%%%
                         for iterPIKI=1:size(luznePIKI,1)
                             testowyRR=nastepnaSZPILA-aktualnaSZPILA-luznePIKI(iterPIKI);
                             fitowyRR = mean(fitTABLE(:,2));
                             tablicaPRZECIEC(iterPIKI)=abs(testowyRR-fitowyRR);
                         end
                         [~,INDEXnowaSZPILA]=min(tablicaPRZECIEC);
                         tester=nastepnaSZPILA-(aktualnaSZPILA+abs(luznePIKI(INDEXnowaSZPILA)));
                         if (tester<upZakresMODA && tester>lowZakresMODA)
                            fills = cat(1,fills,round(nastepnaSZPILA-tester));
                            stopWatch =1; 
                         end
                         if(stopWatch == 0)                       
                            fills = cat(1,fills, round(nastepnaSZPILA-mean(fitTABLE(:,2))));
                            stopWatch =1;
                         end
                     end
                 end      


                 if iter == size(dANNourSTAGE2,1) && Nlength-ANNourSTAGE2(iter+1)>MEANMODA
                    fitTABLE=zeros(4,2);
                     if(iter-1>0)
                        fitTABLE(1,1)=ANNourSTAGE2(iter-1);
                        fitTABLE(1,2)=dANNourSTAGE2(iter-1);
                     end                         
                     if(iter>0)
                        fitTABLE(2,1)=ANNourSTAGE2(iter);
                        fitTABLE(2,2)=dANNourSTAGE2(iter);
                     end                         
                     if(iter+1<=size(dANNourSTAGE2,1))
                        fitTABLE(3,1)=ANNourSTAGE2(iter+1);
                         fitTABLE(3,2)=dANNourSTAGE2(iter+1);
                     end
                     if(iter+2<=size(dANNourSTAGE2,1))
                        fitTABLE(4,1)=ANNourSTAGE2(iter+2);
                        fitTABLE(4,2)=dANNourSTAGE2(iter+2);
                     end
                     fitTABLE(all(fitTABLE==0,2),:)=[];
                     warning('off','all');
                     p = polyfit(fitTABLE(:,1),fitTABLE(:,2),1);
                     warning('on','all'); 
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     aktualnaSZPILA=floor(ANNourSTAGE2(iter));
                     nastepnaSZPILA=Nlength;
                     if aktualnaSZPILA<1
                         aktualnaSZPILA=1;
                     end
                     wycinekSYGNALU=dobrySYGNAL(aktualnaSZPILA:nastepnaSZPILA,1);
                     luznePIKI=findDeccelerations(wycinekSYGNALU,0.8,1,5,40,10000000);              
                     luznePIKI=luznePIKI-1;
                     if(isempty(luznePIKI))
                        if (polyval(p,ANNourSTAGE2(iter))>lowZakresMODA && polyval(p,ANNourSTAGE2(iter))<upZakresMODA)
                            fills = cat(1,fills,round(ANNourSTAGE2(iter+1)+polyval(p,ANNourSTAGE2(iter))));									
                        else
                            fills = cat(1,fills,round(ANNourSTAGE2(iter+1)+MEANMODA));								
                        end
                         stopWatch =1;
                     else                                         
                         tablicaPRZECIEC=zeros(size(luznePIKI,1),1);
                         for iterPIKI=1:size(luznePIKI,1)
                             testowyRR=luznePIKI(iterPIKI);
                             tczs=ANNourSTAGE2(iter);
                             fitowyRR = polyval(p,aktualnaSZPILA);
                             tablicaPRZECIEC(iterPIKI)=abs(testowyRR-fitowyRR);
                         end

                         [~,INDEXnowaSZPILA]=min(tablicaPRZECIEC);
                         tester=abs(luznePIKI(INDEXnowaSZPILA));

                         if (tester<upZakresMODA && tester>lowZakresMODA)
                            fills = cat(1,fills,round(tester+ANNourSTAGE2(iter+1)));
                            stopWatch =1; 
                         end   
                         if(stopWatch == 0)
                            if (polyval(p,ANNourSTAGE2(iter))>lowZakresMODA && polyval(p,ANNourSTAGE2(iter))<upZakresMODA)
                                fills = cat(1,fills,round(ANNourSTAGE2(iter+1)+polyval(p,aktualnaSZPILA)));								
                            else
                                fills = cat(1,fills,round(ANNourSTAGE2(iter+1)+MEANMODA));								
                            end							
                            stopWatch =1;
                         end
                     end                 

                 end


                 if dANNourSTAGE2(iter) > dziuraLimit   
                    
                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     aktualnaSZPILA=floor(ANNourSTAGE2(iter));
                     nastepnaSZPILA=floor(ANNourSTAGE2(iter+1));
                     
                     if aktualnaSZPILA<1
                     elseif nastepnaSZPILA>Nlength
                     else
                         wycinekSYGNALU=dobrySYGNAL(aktualnaSZPILA:nastepnaSZPILA,1);
                         luznePIKI=findDeccelerations(wycinekSYGNALU,lowg-0.3,highg,5,30,10000);
                         wycinekSYGNALUcor=poSkorelowaniuZTemplatemDziecka(aktualnaSZPILA:nastepnaSZPILA,1);
                         
                         luznePIKI2=findDeccelerations(wycinekSYGNALUcor,lowBestPoUsunieciuMatki-0.3,1.0,5,30,10000);   
                         luznePIKI=[luznePIKI; luznePIKI2];
                         luznePIKI=luznePIKI-1;
                          if(isempty(luznePIKI))

%                             if (polyval(p,ANNourSTAGE2(iter))>lowZakresMODA && polyval(p,ANNourSTAGE2(iter))<upZakresMODA)
% 
%                                 fills = cat(1,fills,round(aktualnaSZPILA+polyval(p,ANNourSTAGE2(iter))));
%                                 disp('wypelniam polyvalem bo nie ma luznych pikow');
%                                 disp(num2str(aktualnaSZPILA+polyval(p,ANNourSTAGE2(iter))));
%                                 disp(num2str(polyval(p,ANNourSTAGE2(iter))));
%                             else
% 
%                                 fills = cat(1,fills,round(ANNourSTAGE2(iter)+MEANMODA));		%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX	
%                                 disp('wypelniam MEANMODA bo nie ma luznych pikow a polyval ssie');	
%                                 disp(num2str(ANNourSTAGE2(iter)+MEANMODA));				
%                             end
%                              stopWatch =1;
                         else                                         
                             tablicaPRZECIEC=zeros(size(luznePIKI,1),1);
                             
                               
                            MEANMODAITER = 425;
                            MODAITER = diff(ANNourSTAGE2);
                            if(size(MODAITER,1)>=1)
                                   MODAITER_TABLE=[];
                                   
                               licznik_d=1;
                               iter_d=iter-1; 
                               if iter_d < 1
                                  iter_d=1;
                               end
                                while licznik_d<=4;
                                    if(MODAITER(iter_d) <upZakresMODA && MODAITER(iter_d) >lowZakresMODA )
                                        MODAITER_TABLE(end+1,1) = MODAITER(iter_d);
                                        licznik_d=licznik_d+1;
                                    end                                    
                                    iter_d=iter_d-1;  
                                    if iter_d < 1
                                      break;
                                    end
                                end
                               licznik_g=1;
                               iter_g=iter+1; 
                               if iter_g >= size(MODAITER,1)
                                  iter_g=size(MODAITER,1);
                               end
                                while licznik_g<=4;
                                    if(MODAITER(iter_g) <upZakresMODA && MODAITER(iter_g) >lowZakresMODA )
                                        MODAITER_TABLE(end+1,1) = MODAITER(iter_g);
                                        licznik_g=licznik_g+1;
                                    end   
                                    
                                    iter_g=iter_g+1;  
                                    if iter_g >= size(MODAITER,1)-1
                                      break;
                                    end
                                end                                
                                if(isempty(MODAITER_TABLE))
                                    MODAITER_TABLE = MEANMODAITER;
                                end  
                                %XXXXXXXXXXXXXXXXXXx
                                MEANMODAITER=mean(MODAITER_TABLE);
                            end         
                               
                             for iterPIKI=1:size(luznePIKI,1)
                                 testowyRR=luznePIKI(iterPIKI);
                                 tczs=ANNourSTAGE2(iter);
                                 fitowyRR = MEANMODAITER;
                                 tablicaPRZECIEC(iterPIKI)=abs(testowyRR-fitowyRR);
                             end

                             [~,INDEXnowaSZPILA]=min(tablicaPRZECIEC);
                             tester=abs(luznePIKI(INDEXnowaSZPILA));

                             if (tester<upZakresMODA && tester>lowZakresMODA)
                                fills = cat(1,fills,round(tester+aktualnaSZPILA));
                               % disp('wypelniam luznym pikiem z przodu ');	
                                %disp(num2str(round(tester+aktualnaSZPILA)));				
                                stopWatch =1; 
                             end

                             %%%%%%%%%%%%%%%%
                             for iterPIKI=1:size(luznePIKI,1)
                                 testowyRR=nastepnaSZPILA-aktualnaSZPILA-luznePIKI(iterPIKI);
                                 tczs=nastepnaSZPILA-testowyRR;
                                 fitowyRR = MEANMODAITER;
                                 tablicaPRZECIEC(iterPIKI)=abs(testowyRR-fitowyRR);
                             end
                             [~,INDEXnowaSZPILA]=min(tablicaPRZECIEC);
                             tester=nastepnaSZPILA-(aktualnaSZPILA+abs(luznePIKI(INDEXnowaSZPILA)));
                             if (tester<upZakresMODA && tester>lowZakresMODA)
                                fills = cat(1,fills,round(nastepnaSZPILA-tester));
                               % disp('wypelniam luznym pikiem z tylu ');	
                                %disp(num2str(round(tester+aktualnaSZPILA)));		
                                stopWatch =1; 
                             end
                             if(stopWatch == 0)
                                 if(round(ANNourSTAGE2(iter)+MEANMODAITER)<1)                            
                                 elseif round(ANNourSTAGE2(iter)+MEANMODAITER)>=Nlength
                                 else
                                   % disp('wypelniam fackup');	
                                   % disp(num2str(ANNourSTAGE2(iter)+MEANMODAITER));
                                     fills = cat(1,fills,round(ANNourSTAGE2(iter)+MEANMODAITER));
                                     stopWatch =1;
                                 end
                             end
                         end
                      end
                 end
            end
            ANNourSTAGE2=cat(1,ANNourSTAGE2,fills);
            ANNourSTAGE2=sort(ANNourSTAGE2);
            dANNourSTAGE2=diff(ANNourSTAGE2);
        end
    end
    ANNourSTAGE3=unique(ANNourSTAGE2);
    dANNourSTAGE3 = diff(ANNourSTAGE3);
    
    
    
    for iterFIN=2:size(dANNourSTAGE3)-2
        if dANNourSTAGE3(iterFIN) < 200
            temp=(ANNourSTAGE3(iterFIN-1) + ANNourSTAGE3(iterFIN+2))/2;
            pierwszy= abs(temp-dANNourSTAGE3(iterFIN));
            drugi= abs(temp-dANNourSTAGE3(iterFIN+1));
            if pierwszy < drugi
                 ANNourSTAGE3(iterFIN+1)=999999;
            else
                 ANNourSTAGE3(iterFIN)=999999;
            end
        end
    end
    
    ANNourSTAGE3(ANNourSTAGE3==999999)=[];
    ANNourSTAGE3(ANNourSTAGE3<1)=[];
    ANNourSTAGE3(ANNourSTAGE3>Nlength)=[];
    dANNourSTAGE3=diff(ANNourSTAGE3);    
    
%     ANNourSTAGE3
%     polaczonySzkielet
%     pause;
%       etap ostatni - usun odstaj?ce

     MEANSTD = 30;
     MODA = diff(ANNourSTAGE3);
     if(size(MODA,1)>=1)
         MODA(MODA>upZakresMODA) =[];
         MODA(MODA<lowZakresMODA) =[];
        MEANSTD = 3.2*std(MODA);
     end
     iterSH=1;
     while iterSH < size(dANNourSTAGE3,1)-1
         
         dANNourSTAGE3=diff(ANNourSTAGE3);   
         MEANMODAITER = 425;
         MODAITER = diff(ANNourSTAGE3);
         if (size(MODAITER,1)>=1)
               MODAITER_TABLE=[];

           licznik_d=1;
           iter_d=iterSH-1; 
           if iter_d < 1
              iter_d=1;
           end
            while licznik_d<=5;
                if(MODAITER(iter_d) <upZakresMODA && MODAITER(iter_d) >lowZakresMODA )
                    MODAITER_TABLE(end+1,1) = MODAITER(iter_d);
                    licznik_d=licznik_d+1;
                end                                    
                iter_d=iter_d-1;  
                if iter_d < 1
                  break;
                end
            end
           licznik_g=1;
           iter_g=iterSH+1; 
           if iter_g >= size(MODAITER,1)
              iter_g=size(MODAITER,1);
           end
            while licznik_g<=5;
                if(MODAITER(iter_g) <upZakresMODA && MODAITER(iter_g) >lowZakresMODA )
                    MODAITER_TABLE(end+1,1) = MODAITER(iter_g);
                    licznik_g=licznik_g+1;
                end   

                iter_g=iter_g+1;  
                if iter_g >= size(MODAITER,1)-1
                  break;
                end
            end                                
            if(isempty(MODAITER_TABLE))
                MODAITER_TABLE = MEANMODAITER;
            end  
            %XXXXXXXXXXXXXXXXXXx
            MEANMODAITER=round(mean(MODAITER_TABLE));
        end         
         if  dANNourSTAGE3(iterSH)<=550  && dANNourSTAGE3(iterSH)>=250 && (dANNourSTAGE3(iterSH) >= MEANMODAITER + MEANSTD || dANNourSTAGE3(iterSH) <= MEANMODAITER - MEANSTD  )
                    %MEANMODAITER
                    %MEANSTD
                    dANNourSTAGE3(iterSH);
                   %czy pik w szkielecie jest samotny i do przestawienia?
                    [Lia, Locb]= ismember(round(ANNourSTAGE3(iterSH+1,1)),polaczonySzkielet);
                    if(Locb-1>1 && Locb+1 < size(polaczonySzkielet,1))
                        leftbo = abs(polaczonySzkielet(Locb-1)-polaczonySzkielet(Locb));
                        rightbo = abs(polaczonySzkielet(Locb+1)-polaczonySzkielet(Locb));
                    else
                        leftbo = 9999;
                        rightbo = 9999;
                    end
                    if(Lia && leftbo>550 && rightbo>550 )%XXXXXXXXXXXXXXXXX
                        ANNourSTAGE3(iterSH+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                    elseif(~Lia )
                        ANNourSTAGE3(iterSH+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                       % disp('przestawione bo nie ze szkieletu 1')
                    elseif(~ismember(round(ANNourSTAGE3(iterSH,1)),polaczonySzkielet) )
                        ANNourSTAGE3(iterSH,1)=ANNourSTAGE3(iterSH+1,1)-MEANMODAITER;
                       % disp('przestawione bo nie ze szkieletu 2')
                    end
                    %  disp('przestawione')
         end
               
         if  dANNourSTAGE3(iterSH)>=550 
                %MEANMODAITER
                    ANNourSTAGE3(end+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                      ANNourSTAGE3=sort(ANNourSTAGE3);
                      %disp('dodane');
         end
         if  dANNourSTAGE3(iterSH)<=250
                if  dANNourSTAGE3(iterSH)+dANNourSTAGE3(iterSH+1)>=600                     
                    
                    
                    [Lia, Locb]= ismember(round(ANNourSTAGE3(iterSH+1,1)),polaczonySzkielet);
                    if(Locb-1>1 && Locb+1 < size(polaczonySzkielet,1))
                        leftbo = abs(polaczonySzkielet(Locb-1)-polaczonySzkielet(Locb));
                        rightbo = abs(polaczonySzkielet(Locb+1)-polaczonySzkielet(Locb));
                    else
                        leftbo = 9999;
                        rightbo = 9999;
                    end
                    
                    if(Lia && leftbo>550 && rightbo>550 )%XXXXXXXXXXXXXXXXX
                        ANNourSTAGE3(iterSH+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                    elseif(~Lia )
                        ANNourSTAGE3(iterSH+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                      %  ANNourSTAGE3(iterSH+1,1)
                       % ANNourSTAGE3(iterSH,1)
                       % polaczonySzkielet
                       % pause;
                    %  disp('przestawione bo za krotkie i dziura 1')
                    elseif(~ismember(round(ANNourSTAGE3(iterSH,1)),polaczonySzkielet) )
                        ANNourSTAGE3(iterSH,1)=ANNourSTAGE3(iterSH+1,1)-MEANMODAITER;
                      %disp('przestawione bo za krotkie i dziura 2')
                    end
                else
                    [Lia, Locb]= ismember(round(ANNourSTAGE3(iterSH+1,1)),polaczonySzkielet);
                    if(Locb-1>1 && Locb+1 < size(polaczonySzkielet,1))
                        leftbo = abs(polaczonySzkielet(Locb-1)-polaczonySzkielet(Locb));
                        rightbo = abs(polaczonySzkielet(Locb+1)-polaczonySzkielet(Locb));
                    else
                        leftbo = 9999;
                        rightbo = 9999;
                    end
                    
                    if(Lia && leftbo>550 && rightbo>550 )%XXXXXXXXXXXXXXXXX
                        ANNourSTAGE3(iterSH+1,1)=ANNourSTAGE3(iterSH,1)+MEANMODAITER;
                    elseif(~ismember(round(ANNourSTAGE3(iterSH+1,1)),polaczonySzkielet) )
                        ANNourSTAGE3(iterSH+1,1)=9999999;
                    else
                        ANNourSTAGE3(iterSH,1)=9999999;
                    end
                    ANNourSTAGE3=sort(ANNourSTAGE3);
                    %disp('ODJETE');
                    iterSH=iterSH-1;

                end
         end
         ANNourSTAGE3(ANNourSTAGE3>Nlength)=[];
         iterSH=iterSH+1;
     end
     
     
       
    for iterFIN=2:size(dANNourSTAGE3)-2
        if dANNourSTAGE3(iterFIN) < 250
            temp=(ANNourSTAGE3(iterFIN-1) + ANNourSTAGE3(iterFIN+2))/2;
            pierwszy= abs(temp-dANNourSTAGE3(iterFIN));
            drugi= abs(temp-dANNourSTAGE3(iterFIN+1));
            if pierwszy < drugi
                 ANNourSTAGE3(iterFIN+1)=999999;
            else
                 ANNourSTAGE3(iterFIN)=999999;
            end
        end
    end
    if dANNourSTAGE3(1) < 250
             ANNourSTAGE3(1)=999999;
    end
    if dANNourSTAGE3(size(dANNourSTAGE3,1)-1) < 250
             ANNourSTAGE3(size(dANNourSTAGE3,1))=999999;
    end
    if dANNourSTAGE3(size(dANNourSTAGE3,1)) < 250
             ANNourSTAGE3(size(dANNourSTAGE3,1))=999999;
    end
     
     
    ANNourSTAGE3(ANNourSTAGE3==999999)=[];
    ANNourSTAGE3(ANNourSTAGE3<1)=[];
    ANNourSTAGE3=round(ANNourSTAGE3);
    ANNourSTAGE3(ANNourSTAGE3>Nlength)=[];
    dANNourSTAGE3=diff(ANNourSTAGE3);   
  
    
   %  ---KONIEC--------------------------------------------
   
     %plotting output
% % %      [score1,score2] = readScore(s);
% % %      
% % %     kanal1 = permtable(ktorapermutacja,1);
% % %     kanal2 = permtable(ktorapermutacja,2);
% % %      
% % %      h = figure; 
% % %      subplot(3,1,1);
% % %      set(0,'DefaultAxesColorOrder',[0.6,0.6,0.6])
% % %      %anotacje zewnetrzne
% % %      annfile   = strcat(s(1:end-4),'.fqrs.txt');
% % %      ANN=annotationsFile(annfile, fs);
% % %      ANN(:,1)= ANN(:,1)*1000;
% % %      %nasz szkielet 1.0
% % %      amplitudes3 = ones(1,length(ANNour2)) * 40;
% % %      amplitudes3 = amplitudes3';
% % %      ANNestSZKIELET = [ANNour2,amplitudes3];
% % % %     nasz szkielet 1.0 po?redni
% % %      amplitudes = ones(1,length(ANNourREF)) * 55;
% % %      amplitudes = amplitudes';
% % %      ANNest = [ANNourREF,amplitudes];
% % % %     nasz szkielet 2.0 est
% % %      amplitudes6 = ones(1,length(ANNourREF20)) * 55;
% % %      amplitudes6 = amplitudes6';
% % %      ANNest_20 = [ANNourREF20,amplitudes6];
% % %      %nasz szkielet polaczony
% % %      amplitudes8 = ones(1,length(polaczonySzkielet)) * 55;
% % %      amplitudes8 = amplitudes8';
% % %      ANNpolaczonySzkielet = [polaczonySzkielet,amplitudes8];
% % %      
% % %      
% % %      %nasze  ostateczne oznaczenia
% % %      amplitudes9 = ones(1,length(ANNourSTAGE3)) * 60;
% % %      amplitudes9 = amplitudes9';
% % %      ANNostateczneUZU = [ANNourSTAGE3,amplitudes9];
% % %      %nasz szkielet 2.0
% % %      amplitudes7 = ones(1,length(ANNour_szkielt20)) * 40;
% % %      amplitudes7 = amplitudes7';
% % %      ANNestSZKIELET_20 = [ANNour_szkielt20,amplitudes7];
% % %      
% % %      hold on;
% % %      plot(dobrySYGNAL(:,1));
% % %      stem(ANN(:,1),ANN(:,2),'or');
% % %      if(isempty(ANNest))
% % %      else
% % %             stem(ANNestSZKIELET(:,1),ANNestSZKIELET(:,2),'ok');
% % %      end     
% % %      if(isempty(ANNostateczneUZU))
% % %      else
% % %             stem(ANNostateczneUZU(:,1),ANNostateczneUZU(:,2),'ob');
% % %          
% % %      end
% % %      axis([0,60000,-80,80]);
% % %      legend({'ECG', 'original ANN', 'our ANN est','ostateczne ANN'},'Location','SouthWest');
% % %      ccc = [ 'WYNIK 1(HR):', num2str(score1), ' WYNIK 2(RR):', num2str(score2), ' granice: dolna:' num2str(lowg) ', gorna:' num2str(highg) ', dlugosc_d:' num2str(lowLengthBest) ', dlugosc_g:' num2str(highLengthBest) ', zero:',num2str(zeroBest),', kanaly:' num2str(kanal1) ', ' num2str(kanal2) ', ktory wziety' num2str(whichchannel) ];
% % %      title(ccc);
% % %      ylabel('ECG');
% % %      xlabel('time[samples]');
% % %      
% % %      
% % %      subplot(3,1,2);
% % %       hold on;
% % %      
% % %       plot(odjetaMatka(:,1));
% % %        plot(signalZeroMother(:,1),'or','MarkerSize',2);
% % %       stem(ANN(:,1),ANN(:,2),'or');
% % %       if(isempty(ANNpolaczonySzkielet))
% % %       else
% % %             stem(ANNpolaczonySzkielet(:,1),ANNpolaczonySzkielet(:,2),'ok');
% % %       end    
% % %      if(isempty(ANNostateczneUZU))
% % %      else
% % %             stem(ANNostateczneUZU(:,1),ANNostateczneUZU(:,2),'ob');
% % %          
% % %      end
% % %       axis([0,60000,-80,80]);
% % %       legend({'ECG po odjeciu MATKI','korelacja dziecko', 'original ANN', 'our ANN2 est', 'ostateczne ANN'},'Location','SouthWest');
% % %      ccc = [ 'WYNIK 1:', num2str(score1), ' WYNIK 2:', num2str(score2),' granice: dolna:' num2str(lowg) ', gorna:' num2str(highg) ', dlugosc_d:' num2str(lowLengthBest) ', dlugosc_g:' num2str(highLengthBest) ', zero:',num2str(zeroBest) ', kanaly:' num2str(kanal1) ', ' num2str(kanal2) ', ktory wziety' num2str(whichchannel) ];
% % %      title(ccc);
% % %       ylabel('ECG');
% % %       xlabel('time[samples]');
% % %      
% % %       
% % %      subplot(3,1,3);
% % %      
% % %      
% % %      rryzICH=diff(ANN(:,1));
% % %      positionRRyzICH=ANN(1:size(ANN,1)-1,1);
% % %      rryNaszeOST=diff(ANNostateczneUZU(:,1));
% % %      positionrryNaszeOST=ANNostateczneUZU(1:size(ANNostateczneUZU,1)-1,1);
% % %      
% % %      plot(positionRRyzICH,rryzICH,'r');
% % %      hold on;
% % %      plot(positionrryNaszeOST,rryNaszeOST,'ok');
% % % 
% % %         
% % %      outimg1 = strcat(s(1:end-4),'.dat.png');
% % %      set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
% % %      print(h,'-r300','-dpng',outimg1);
% % %   %  close(h);    
% % % 
% % %      h2 = figure(2); 
% % % 
% % %      plot(templateQT(:,1));
% % %      hold on;
% % %      annQT=[positionOfQ;positionOfT];
% % %      amplitudesQT = ones(1,length(annQT)) * 0.1;
% % %      amplitudesQT = amplitudesQT';
% % %      ANNestQT= [annQT,amplitudesQT];
% % %       if(isempty(ANNestQT))
% % %       else
% % %             stem(ANNestQT(:,1),ANNestQT(:,2),'or');
% % %       end
% % %       
% % %      ccc = [ 'QT interval:' num2str(ourQT) ];
% % %      title(ccc);
% % %       
% % %      outimg2 = strcat(s(1:end-4),'.QT.png');
% % %      set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 9])
% % %      print(h2,'-r300','-dpng',outimg2);
    % close(h2);    

    %---KONIEC--------------------------------------------
     fetal_QRSAnn_est    = ANNourSTAGE3;
     QT_Interval         = ourQT;
     
     
end

function [y]=wypelniaczPustek(tm,ECG)
    N       = size(ECG,2);     % number of abdominal channels
    y=ECG;
    for iterPoKanalach=1:N
        %wez 1 kanal ECG
        x = ECG(:,iterPoKanalach);
        %wyszukaj indeksy pozycji z NaN;
        pustki = find(isnan(x) == 1);        
        %Jesli sa puste
        if isempty(pustki)~=1
            dggg=zeros(2,1);
            iloscPustek=1;
            if size(pustki,1) >= 2
                dg=1;
                for iter0=1:(size(pustki,1)-1)
                    if pustki(iter0+1) - pustki(iter0) ~= 1
                        dggg(1,iloscPustek)=pustki(dg);
                        dggg(2,iloscPustek)=pustki(iter0);
                        dg=iter0+1;
                        iloscPustek=iloscPustek+1;
                    end
                end
            end
            if size(pustki,1) >= 2
                        dggg(1,iloscPustek)=pustki(dg);
                        dggg(2,iloscPustek)=pustki(size(pustki,1));
            end
            if size(pustki,1) == 1
                dggg=[pustki(1) pustki(1)] ;
            end
            
            % a teraz parabolowanie po wszystkich dziurach
            for iter201=1:size(dggg,2);
                dziurkacz = dggg(:,iter201);
                dziurapom = zeros(10,2);
                for itertemp=1:5
                    if (dziurkacz(1)-6+itertemp) > 0
                        dziurapom(itertemp,1) = tm(dziurkacz(1)-6+itertemp);
                        dziurapom(itertemp,2) = ECG(dziurkacz(1)-6+itertemp,iterPoKanalach);
                    end
                end
                for itertemp=6:10
                    if (dziurkacz(2)+itertemp-5) <= size(tm,1)
                        dziurapom(itertemp,1) = tm(dziurkacz(2)+itertemp-5);
                        dziurapom(itertemp,2) = ECG(dziurkacz(2)+itertemp-5,iterPoKanalach);
                    end
                end
                dziurapom(all(dziurapom==0,2),:)=[];
                warning('off','all');
                p = polyfit(dziurapom(:,1),dziurapom(:,2),2);
                warning('on','all');
                for itertemp=dziurkacz(1):dziurkacz(2)
                    ECG(itertemp,iterPoKanalach)=p(1)*tm(itertemp)^2 + p(2)*tm(itertemp)+ p(3);
                    y(itertemp,iterPoKanalach)=p(1)*tm(itertemp)^2 + p(2)*tm(itertemp)+ p(3);
                    
                end
            end
        end
    end
end


function [skeleton]=skeletonJoin(ske_base,ske_new,lowZakresMODA)

    starySzkielet = ske_base;
    nowySzkielet = ske_new;
    if(size(nowySzkielet,1)>1)
        for iterStary=1:size(starySzkielet,1)
            for iterNowy=1:size(nowySzkielet,1)
                if( abs(starySzkielet(iterStary,1) - nowySzkielet(iterNowy,1)) < lowZakresMODA )
                    nowySzkielet(iterNowy,1) = 999999;
                end
            end
        end
    end
    nowySzkielet(nowySzkielet==999999)=[];
    
    skeleton=[starySzkielet;nowySzkielet];
    skeleton = sort(skeleton);

end


function [cleanSkeleton] = annCleaning(tablicaWykryc, upZakresMODA,lowZakresMODA,oknoDoboruPiku,skeletonLimit)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % KONIEC - MAMY NAJLEPIEJ ZNALEZIONE PIKI
    % teraz - dobor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %wylicz najczesciej pojawiajaca sie wartosc miedzy 350 a 500
    
    MEANMODA = 425;
    %MEANSTD = 30;
    MODA = diff(tablicaWykryc);
    if(size(MODA,1)>=1)
        MODA(MODA>upZakresMODA) =[];
        MODA(MODA<lowZakresMODA) =[];
        MEANMODA=mean(MODA);
        %MODA(MODA>MEANMODA+skeletonLimit) =[];
       % MODA(MODA<MEANMODA-skeletonLimit) =[];
       % MEANMODAdoSTD=mean(MODA)
        %MEANSTD = 0.8*std(MODA)
    end        
    %etap pierwszy - znajdz dobre pobudzenia
    ANNour = tablicaWykryc;
    dANNour = diff(ANNour);
    
    ANNour(1)=999999;
    ANNour(size(ANNour,1))=999999;
    if(size(dANNour,1)>1)
        for iter=2:size(dANNour,1)
            if dANNour(iter) < MEANMODA+skeletonLimit && dANNour(iter) > MEANMODA-skeletonLimit &&...
               dANNour(iter-1) < MEANMODA+skeletonLimit && dANNour(iter-1) > MEANMODA-skeletonLimit
            else
                ANNour(iter)=999999;
            end
        end
    end
    ANNour2=ANNour;
    %etap drugi - znajdz dobre na granicy    
    for iter=1:size(ANNour,1)-1
        if ANNour(iter) == 999999 && ANNour(iter+1) ~= 999999
%             %XXXXX dlugosc
             if(iter+1>size(dANNour,1))
                   referencyjnyRR=MEANMODA;
             else
                 referencyjnyRR = dANNour(iter+1);
             end
             referencyjnyPIK = ANNour(iter+1);
             closematch=[];                
             for iterwew=1:size(tablicaWykryc,1)
                 temp= referencyjnyPIK-tablicaWykryc(iterwew);
                if temp <=referencyjnyRR+oknoDoboruPiku   && temp >=referencyjnyRR-oknoDoboruPiku 
                     closematch(end+1,1) = abs(temp-referencyjnyRR);
                     closematch(end,2) = iterwew;
                end
            end
            if isempty(closematch)
            else
                 [~, iii]=min(closematch(:,1));
                 tablicaWykryc(closematch(iii,2));
                 ANNour2 = cat(1,ANNour2,tablicaWykryc(closematch(iii,2)));
            end
        end
    end   
    for iter=1:size(ANNour,1)-1
        if ANNour(iter+1) == 999999 && ANNour(iter) ~= 999999
            %XXXXX dlugosc
             if(iter-1>0)
                 referencyjnyRR = dANNour(iter-1);
             else
                   referencyjnyRR=MEANMODA;
             end
            referencyjnyPIK = ANNour(iter);
            closematch=[];                
            for iterwew=1:size(tablicaWykryc,1)
                temp= tablicaWykryc(iterwew)-referencyjnyPIK;
                if temp <=referencyjnyRR+oknoDoboruPiku   && temp >=referencyjnyRR-oknoDoboruPiku  
                     closematch(end+1,1) = abs(temp-referencyjnyRR);
                     closematch(end,2) = iterwew;
                end
            end
            if isempty(closematch)
            else
                 [~, iii]=min(closematch(:,1));
                  ANNour2 = cat(1,ANNour2,tablicaWykryc(closematch(iii,2)));
            end 
        end
    end        
    ANNour2(ANNour2==999999)=[];
    ANNour2 = sort(ANNour2);
    dANNour2=diff(ANNour2);
    
    %etap trzeci - usun podwojne na granicy
    for iterSH=2:size(dANNour2)-2
        if dANNour2(iterSH) < 2*(upZakresMODA - lowZakresMODA)
            temp=(ANNour2(iterSH-1) + ANNour2(iterSH+2))/2;
            pierwszy= abs(temp-dANNour2(iterSH));
            drugi= abs(temp-dANNour2(iterSH+1));
            if pierwszy < drugi
                 ANNour2(iterSH+1)=999999;
            else
                 ANNour2(iterSH)=999999;
            end
        end
    end
    
    ANNour2(ANNour2==999999)=[];
    

    cleanSkeleton = ANNour2;
end

function [firstBest, secondBest] = brutforsPorownanie(tablicaWykryc, upZakresMODA,lowZakresMODA,rangeSD,skeletonLimit)

    tabliceSzkieletow3=zeros(size(tablicaWykryc(:,1),1),1);
    for iter=1:size(tablicaWykryc(:,1))
        A = cell2mat(tablicaWykryc(iter,3));
        MODA = diff(A);
        if(size(MODA,1)>=1)
            MODA(MODA>upZakresMODA) =[];
            MODA(MODA<lowZakresMODA) =[];
        end        
        if(isempty(MODA))
            tabliceSzkieletow3(iter)=0;
        else
            tabliceSzkieletow3(iter)=size(MODA,1);
        end
        MEANMODA = mean(MODA);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ANNour = A;
        dANNour =diff(A);
        if(size(dANNour,1)>1)
            for iterS=2:size(dANNour,1)
                %XXXX zacina sie przy dlugosci = 1
                if dANNour(iterS) < MEANMODA+skeletonLimit && dANNour(iterS) > MEANMODA-skeletonLimit &&...
                   dANNour(iterS-1) < MEANMODA+skeletonLimit && dANNour(iterS-1) > MEANMODA-skeletonLimit
                else
                    ANNour(iterS)=999999;
                end
            end
        end
        ANNour(ANNour==999999)=[];
        ANNour = sort(ANNour);
        if(isempty(ANNour))
            tabliceSzkieletow3(iter)=0;
        else
            tabliceSzkieletow3(iter)=size(ANNour,1);
        end
        %----warnuek na SD-------------------------------------------------
        if std(dANNour) > rangeSD
            tabliceSzkieletow3(iter) = NaN;
        end         
    end    
    [~,I] = max(tabliceSzkieletow3);    
    firstBest = cell2mat(tablicaWykryc(I(1,1),1));
    secondBest = cell2mat(tablicaWykryc(I(1,1),2));
end


function [y]=notchFilter(ECG,frequency)
    fs = 1000;             
    f0 = frequency;               
    fn = fs/2;            
    freqRatio = f0/fn;      %ratio of notch freq. to Nyquist freq.

    notchWidth = 0.15;       %width of the notch

    zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%Compute poles
    poles = (1-notchWidth) * zeros;

   % figure;
   % zplane(zeros.', poles.');

    b = poly( zeros ); % Get moving average filter coefficients
    a = poly( poles ); % Get autoregressive filter coefficients

   % figure;
    %freqz(b,a,32000,fs)

    y = filter(b,a,ECG);
end


function [y]=movingMedian(ECG,window)
    y=ECG;
    if window==1
        return
    end
    N = size(ECG,2);
    for iter1=1:N;
        %[[1]] wez 1 kanal ECG
              x = ECG(:,iter1);
        %[[2]] wez utworz wektor wynikowy
              ECGtemp=zeros(size(x));
        %[[3]] uzupelnij wektor wejsciowy o poczatek i koniec zerami
              x=[zeros(floor(window*.5),1);x;zeros(floor(window*.5),1)];
              window1=window-1;
        %[[5]] wylicz mediane
              for ii=1:size(ECGtemp,1);
                 ECGtemp(ii,:)=median(x(ii+(0:window1),:));
              end
        %[[6]] wstaw do wyniku
              y(:,iter1) = ECG(:,iter1) - ECGtemp;
    end
end

function [y]=pearsonCovariance(ECG,template,przesuniecie)
    x = ECG;
    window=size(template,1);    
    %[[2]] wez utworz wektor wynikowy
        y=zeros(size(ECG));
    %[[3]] uzupelnij wektor wejsciowy o poczatek i koniec zerami
        x=[zeros(przesuniecie,1);x;zeros(window-przesuniecie,1)];
        %x=[x;zeros(window,1)];
        window1=window-1;
    %[[5]] wylicz kowariancje
        for ii=1:size(y,1);
                 C=cov(x(ii+(0:window1),1),template);
                 y(ii,1)=C(2);%/(std(x(ii+(0:window1),1))*std(template));
                 %y(ii,1) = y(ii,1);%*y(ii,1);
        end
    
end
    
function [annotations] = findAnnotationsByTreshold(sygnal, treshold,minimal_distance_between_ann)
    sortowane = sort(sygnal);
   % sortowane
    %disp('size(sygnal)');
    %size(sygnal)
   % disp('round(size(sortowane,1)*treshold)');
   % round(size(sortowane,1)*treshold)
    prog = sortowane(round(size(sortowane,1)*treshold));
    %sortowane(round(size(sortowane,1)*treshold)-1)
    % XXXnie omijamy pierwszego i ostatniego
    tableOfOneRun=[];
    annotations=[];
    %odnajdywannie ci?gów z maximami
    iterator = 1;
    while ( iterator  < size(sygnal,1) )
        if( sygnal(iterator,1) >= prog)
            tableOfOneRun(end+1,1) = iterator;     
            tableOfOneRun(end,2) = sygnal(iterator,1);
        else
            if( ~isempty(tableOfOneRun) )
                [C,I]= max(tableOfOneRun(:,2));
                annotations(end+1,1)=tableOfOneRun(I,1);
                iterator=iterator+minimal_distance_between_ann;
            end
            tableOfOneRun=[];
        end 
        iterator=iterator+1;
    end
end


function [templateChild,ourQT,positionOfQ,positionOfT]=estimateQT(ECG, szkielet)

    zakresLewy = -40;
    zakresPrawy = +390;
    % wyznaczenie template'u dziecka
    tableOfChildrens = zeros( zakresPrawy-zakresLewy+1, size(szkielet,1));
    templateChild = zeros( zakresPrawy-zakresLewy+1, 1);
    if(size(szkielet,1)>1 && size(szkielet,2)>0)
        for iterM=1:size(szkielet,1)
            %XXXXXmozliwe dzielenie przez 0
            if( szkielet(iterM,1)+zakresLewy >=1 && (szkielet(iterM,1)+zakresPrawy < size(ECG,1)-1))
                tableOfChildrens(:,iterM) = ECG( szkielet(iterM,1)+zakresLewy:szkielet(iterM,1)+zakresPrawy,1)/abs(ECG(szkielet(iterM,1),1));
            end
        end
    end
    templateChild = mean(tableOfChildrens,2);
    childPeakIndex = -zakresLewy+1;
    
    
    %[~,positionOfQ_min] = min(templateChild(1:childPeakIndex,1));
    [~,positionOfQ] = max(templateChild(5:25,1));
    positionOfQ = positionOfQ + 4;
    
    [~,positionOfT] = min(templateChild(positionOfQ+150:positionOfQ+300,1));
    positionOfT=positionOfT+positionOfQ+150;
    ourQT=positionOfT-positionOfQ;
    
end

function [korelacjaDziecko]=templateOfChild(ECG, szkielet)
    zakresLewy = -25;
    zakresPrawy = +35;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wyznaczenie template'u dziecka
    tableOfChildrens = zeros( zakresPrawy-zakresLewy+1, size(szkielet,1));
    templateChild = zeros( zakresPrawy-zakresLewy+1, 1);
    for iterM=1:size(szkielet,1)
        %XXXXXmozliwe dzielenie przez 0
        if( szkielet(iterM,1)+zakresLewy >=1 && (szkielet(iterM,1)+zakresPrawy < size(ECG,1)-1))
            tableOfChildrens(:,iterM) = ECG( szkielet(iterM,1)+zakresLewy:szkielet(iterM,1)+zakresPrawy,1)/abs(ECG(szkielet(iterM,1),1));
        end
    end
    templateChild = mean(tableOfChildrens,2);
    
    childPeakIndex = -zakresLewy+1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
    korelacjaDziecko = pearsonCovariance(ECG,templateChild(childPeakIndex+zakresLewy:childPeakIndex+zakresPrawy,1),-zakresLewy);
    %%%XXXXXXXXXXgranica do kowariancji na sztywno
   
end

function [MotherEcgAnnotations, ecgMother,korelacjaMatka,anotacjePrzedMatka]=usuwaczMatki(ECG)
    % zakresy template'u
    %zakresLewy = -80;
    %zakresPrawy = +100;
    zakresLewy = -140;
    zakresPrawy = +280;
    % zakres do kowariancji
    zakresLewyCovariancja = -20;
    zakresPrawyCovariancja = 70;
    
    prog_templateu995 = 0.998;
    prog_templateu990 = 0.995;
    prog_templateu985 = 0.990;
    prog_templateu980 = 0.985;
    minimalna_odl_miedzy_matkami = 200;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    anotacjePrzedMatka = findAnnotationsByTreshold(abs(ECG), prog_templateu995,minimalna_odl_miedzy_matkami); 

    anotacjePrzedMatkaNEW = findAnnotationsByTreshold(abs(ECG), prog_templateu990,minimalna_odl_miedzy_matkami);
    if(min(diff(anotacjePrzedMatkaNEW))<400)
        %disp('poziom 995');    
    else
        anotacjePrzedMatka = anotacjePrzedMatkaNEW;
        anotacjePrzedMatkaNEW=findAnnotationsByTreshold(abs(ECG), prog_templateu985,minimalna_odl_miedzy_matkami); 
        if(min(diff(anotacjePrzedMatkaNEW))<400)
           % disp('poziom 990');    
        else
            anotacjePrzedMatka = anotacjePrzedMatkaNEW;
            anotacjePrzedMatkaNEW=findAnnotationsByTreshold(abs(ECG), prog_templateu980,minimalna_odl_miedzy_matkami);  
            if(min(diff(anotacjePrzedMatkaNEW))<400)
              % disp('poziom 985');   
            else
                anotacjePrzedMatka = anotacjePrzedMatkaNEW;
               % disp('poziom 980');   
            end
        end
    end
    size(anotacjePrzedMatka);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wyznaczenie template'u matki
    tableOfMothers = zeros( zakresPrawy-zakresLewy+1, size(anotacjePrzedMatka,1));
    templateMother = zeros( zakresPrawy-zakresLewy+1, 1);
    size(templateMother,1);
    for iterM=1:size(anotacjePrzedMatka,1)
        %XXXXXmozliwe dzielenie przez 0
        if( anotacjePrzedMatka(iterM,1)+zakresLewy >=1 && (anotacjePrzedMatka(iterM,1)+zakresPrawy < size(ECG,1)-1))
            tableOfMothers(:,iterM) = ECG( anotacjePrzedMatka(iterM,1)+zakresLewy:anotacjePrzedMatka(iterM,1)+zakresPrawy,1)/abs(ECG(anotacjePrzedMatka(iterM,1),1));
        end
    end
    templateMother = mean(tableOfMothers,2);
    motherPeakIndex = -zakresLewy+1;
    %mother Peak Index brany z granic, nie z wykrycia
    %[C,motherPeakIndex]= max(templateMother(:,1));
    % motherPeakIndex
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    korelacjaMatka = pearsonCovariance(ECG,templateMother(motherPeakIndex+zakresLewyCovariancja:motherPeakIndex+zakresPrawyCovariancja,1),-zakresLewyCovariancja);
    %%%XXXXXXXXXXgranica do kowariancji na sztywno
    %TUTAJ
    MotherEcgAnnotations=[];
    MotherEcgAnnotations=findAnnotationsByTreshold(korelacjaMatka, 0.98, minimalna_odl_miedzy_matkami);
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TWORZENIE EKG MATKI Z TEMPLATE'U
    
    tableOfMultipliers = zeros( size(MotherEcgAnnotations,1),1);
    for iterM=1:size(MotherEcgAnnotations,1)
        %XXXX?rednia z wi?kszego zakresu
        %tableOfMultipliers(iterM,1) = mean(abs(ECG( MotherEcgAnnotations(iterM,1)-motherPeakIndex:MotherEcgAnnotations(iterM,1)-motherPeakIndex+size(templateMother,1),1)));
        %tableOfMultipliers(iterM,1) = tableOfMultipliers(iterM,1) / mean(abs(templateMother(:,1)));
        %XXXX?rednia z mniejszego zakresu
        if( MotherEcgAnnotations(iterM,1)+zakresLewyCovariancja >=1 && (MotherEcgAnnotations(iterM,1)+zakresPrawyCovariancja < size(ECG,1)-1))
            tableOfMultipliers(iterM,1) = mean(abs(ECG( MotherEcgAnnotations(iterM,1)+zakresLewyCovariancja:MotherEcgAnnotations(iterM,1)+zakresPrawyCovariancja,1)));
            tableOfMultipliers(iterM,1) = tableOfMultipliers(iterM,1) / mean(abs(templateMother(motherPeakIndex+zakresLewyCovariancja:motherPeakIndex+zakresPrawyCovariancja,1)));
        else
            tableOfMultipliers(iterM,1) = 0;
        end
    end
    ecgMother = zeros(size(ECG,1),1);
    %size(tableOfMultipliers)
    %size(ecgMother( MotherEcgAnnotations(1,1)-motherPeakIndex:MotherEcgAnnotations(1,1)-motherPeakIndex+size(templateMother,1)-1,1))
    for iterM=1:size(MotherEcgAnnotations,1)
        if( MotherEcgAnnotations(iterM,1)-motherPeakIndex+1 >=1 && (MotherEcgAnnotations(iterM,1)-motherPeakIndex+size(templateMother,1) < size(ECG,1)-1))
            ecgMother( MotherEcgAnnotations(iterM,1)-motherPeakIndex+1:MotherEcgAnnotations(iterM,1)-motherPeakIndex+size(templateMother,1),1) = templateMother * tableOfMultipliers(iterM,1);%ampMeanMother;%
        end
    end
end

function [signal] = zeroMother(motherAnn, ECG, zakresLewy, zakresPrawy)
    signal = ECG;
    for iterM=1:size(motherAnn,1)
        %XXXXXmozliwe dzielenie przez 0
        if( motherAnn(iterM,1)+zakresLewy >=1 && (motherAnn(iterM,1)+zakresPrawy < size(ECG,1)-1))

             signal( motherAnn(iterM,1)+zakresLewy+1:motherAnn(iterM,1)+zakresPrawy,1) =0;%ampMeanMother;%
     
        end
    end
end
function annotations = findDeccelerations(ECGstage2,lowBorderDist,upBorderDist,lowLength,upLength,czyZero)
    N=size(ECGstage2,1);
    sortowane = sort(abs(ECGstage2));
    prog1 = sortowane(round(size(sortowane,1)*lowBorderDist));
    if(upBorderDist == 1)
        prog2 = max(sortowane) + 1.0;
    else
        prog2 = sortowane(round(size(sortowane,1)*upBorderDist));
    end
    prog3 = prog1/50;
    
    ECGdiffs=diff(ECGstage2);
    SplitPositions=[];
    
    iter4 = 1;
    while ( iter4 < size(ECGdiffs,1)-1 )
        if ECGdiffs(iter4,1) <= 0 && ECGdiffs(iter4+1,1) <= 0 
            iter4=iter4+1;
            
        elseif ECGdiffs(iter4,1) <= 0 && ECGdiffs(iter4+1,1) <= prog3 
            if ECGdiffs(iter4+2,1) >0
                SplitPositions(end+1,1)=iter4+1;
                iter4=iter4+1;
            else
                iter4=iter4+2;
            end
        else             
            SplitPositions(end+1,1)=iter4;
            iter4=iter4+1;
        end   
    end
    
    piki = [];
    if isempty(SplitPositions)~=1
        dlugoscipochodow = diff(SplitPositions);
        for iter5=1:size(dlugoscipochodow,1)-1
            if lowLength <= dlugoscipochodow(iter5,1) &&...
               upLength  >= dlugoscipochodow(iter5,1) &&...
               prog1 <= ECGstage2(SplitPositions(iter5),1) &&...
               ECGstage2(SplitPositions(iter5+1)-1,1) <= czyZero &&...
               prog2 >= ECGstage2(SplitPositions(iter5),1)
                     MAXneighbour=0;
                     if(SplitPositions(iter5)-20>=1 &&SplitPositions(iter5)+20 <= N)
                           for iterRR=SplitPositions(iter5)-20:SplitPositions(iter5)-5
                             if  abs(ECGstage2(iterRR))>MAXneighbour
                                 MAXneighbour= max(abs(ECGstage2(iterRR)));
                             end
                           end
                           for iterRR=SplitPositions(iter5)+5:SplitPositions(iter5)+20
                             if  abs(ECGstage2(iterRR))>MAXneighbour
                                 MAXneighbour= max(abs(ECGstage2(iterRR)));
                             end
                           end
                     end
                     if(MAXneighbour>prog2)
                     else
                        piki(end+1,1) = SplitPositions(iter5);
                     end
            end  
        end
    end
    annotations=piki;
end

function ANN = annotationsFile(filename,fsample)
    %wczytaj anotacj z pliku
    Amatrix = dlmread(filename);
    %przelicz na czas
    Amatrix = Amatrix / fsample;
    %generacja amplitud dla czasu = 50;
    amplitudes = ones(1,length(Amatrix)) * 70;
    %tablica wierszowa na kolumne
    amplitudes = amplitudes';
    %tablica 2xN czas:amplituda dla anotacji - na wyjscie
    ANN = [Amatrix,amplitudes];
end

function [score1,score2] = readScore(s)
    x = str2num(s(2:end-4));
    %wczytaj anotacj z pliku
    SCORES = dlmread('PCscores.txt');
    
    ind = find(SCORES(:,1)==x);
    score1 = SCORES(ind(1,1),3);
    score2 = SCORES(ind(1,1),4);
end