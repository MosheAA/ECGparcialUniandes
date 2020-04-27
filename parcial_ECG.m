%%% Parcial de ECG - -
% Moshé Alonso Amarillo 
% Universidad Nacional de Colombia - Universidad de los Andes
% 2020
%% Importación de datos 

files = ["001N1_ECG.mat", "002N7_ECG.mat", "004N8_ECG.mat",...  % Lista de archivos 
"005N1_ECG(1).mat","006N1_ECG.mat","007N1_ECG.mat","009N7_ECG.mat"];

v_ECGSig_Cell {1,length(files)}     = [];                       % Preasignación de células 
v_HypCode_Cell {1,length(files)}    = [];                       % Preasignación de células 
for i = 1:length(files)
    load(files(i))    
    v_ECGSig_Cell{1,i} = v_ECGSig;
    v_HypCode_Cell{1,i}= v_HypCode;
    clear v_ECGSig 
    clear v_HypCode
end                                      % Cargar archivos
clear files i 

%% Descripción de la base de datos: 
% s_FsHz          = Frecuencia de muestreo del ECG
% v_ECGSig_Cell   = Señal de ECG                                           
% v_HypCode_Cell  = Hipnograma 
% v_HypTime       = Tiempo de muestreo del hipnograma (Fs = 0.03Hz)
% v_HypCodeLabels = Etiquetas de cada estado en el hipnograma

%%% Observaciones: 
% 1. Las señales no tienen componentes de ruido de línea.
% 2. Las señales tienen distintas duraciones (Al rededor de 8 horas). 

%% Obtención de picos -R- Algoritmo Pan-Tompkins 
% 1. Derivación
for i = 1:length(v_ECGSig_Cell)
    v_ECGFiltDiff{1,i} = zeros (1, numel(v_ECGSig_Cell{1,i}));                      
end                              % Preasignación de células 
for i = 1:length(v_ECGSig_Cell)
    v_ECGFiltDiff_temp          = (zeros (1, numel(v_ECGSig_Cell{1,i})));% Vector temporal
    v_ECGFiltDiff_temp(2:end)   = diff(v_ECGSig_Cell{1,i});     % Derivada ( x'[n]=x[n]-x[n-1])
    v_ECGFiltDiff_temp(1)       = v_ECGFiltDiff_temp(2);        % Se pierde un punto x[1]=x[2]
    v_ECGFiltDiff {1,i}         = v_ECGFiltDiff_temp(1:end);    % Incluir vector temporal en la célula
    clear v_ECGFiltDiff_temp
end                              % 1° Derivada de la señal
% 2. Elevar al cuadrado
for i = 1:length(v_ECGSig_Cell)
    v_ECGFiltDiffsqred{1,i} = zeros (1, numel(v_ECGSig_Cell{1,i}));                      
end                              % Preasignación de células 
for i = 1:length(v_ECGSig_Cell)
    v_ECGFiltDiffsqred_temp  = (v_ECGFiltDiff{1,i}).^2;         % Elevar al cuadrado
    v_ECGFiltDiffsqred {1,i} = v_ECGFiltDiffsqred_temp(1:end);  % Incluir vector temporal en la célula
    clear v_ECGFiltDiffsqred_temp
end                              % Elevar al cuadrado la señal
% 3. Integral móvil 
s_QRSdursec    = 0.03;                                          % Tiempo promedio de QRS
s_AvePointsNum = round(s_QRSdursec*s_FsHz);                     % Ventana de integración en muestras
v_ecginteg {1,length(v_ECGSig_Cell)} = [];
for i = 1:length(v_ECGSig_Cell)
    v_ecginteg{1,i} = zeros (1, numel(v_ECGSig_Cell{1,i}));                     
end                              % Preasignación de células 
for i = 1:length(v_ECGSig_Cell)
    v_ecginteg_temp          = (zeros (1, numel(v_ECGSig_Cell{1,i}))); % Vector temporal integral 
    v_ECGFiltDiffsqred_temp  =  v_ECGFiltDiffsqred {1,i};              % Vector temporal derivada 
    for s_ind = 1:numel(v_ecginteg_temp)
        s_firstind = s_ind -(s_AvePointsNum - 1);
        if s_firstind < 1
            s_firstind = 1;    
        end 
        v_ecginteg_temp(s_ind) = mean(v_ECGFiltDiffsqred_temp(s_firstind:s_ind));
    end                            % Integral 
    v_ecginteg{1,i} = v_ecginteg_temp;                                      % Incluir vector temporal en la célula
    clear v_ecginteg_temp v_ECGFiltDiffsqred_temp s_ind s_firstind
end                              % Integral móvil de la señal
% 4. Detección de -R-
v_peaksind_Cell{1,length(v_ECGSig_Cell)} = [];                  % Preasignación de células
v_peaks_Cell{1,length(v_ECGSig_Cell)}    = [];           
     % Parámetros de umbral para detección de -R-
s_QRSInterDurSec = 0.2;                                         % T. max. entre picos = 0.2 segundos
s_mindursam      = round(s_QRSInterDurSec*s_FsHz);              % Tiempo -> Muestras
min_amp          = 1;                                           % Amplitud mínima = +1 SD
max_amp          = 10;                                           % Amplitud máxima = +8 SD
for i = 1:length(v_ECGSig_Cell)
    v_ecginteg_temp  = v_ecginteg{1,i}; 
   
    v_peaks          = findpeaks (v_ecginteg_temp);             % Hallar picos preliminares 
    s_peaksMean      = mean(v_peaks);                           % Media de los picos detectados
    s_peaksSD        = std(v_peaks);                            % Desviación estándar de los picos 
    s_minThreshold   = s_peaksMean +min_amp*s_peaksSD;          % Amplitud mínima = +1 SD
    s_maxThreshold   = s_peaksMean +max_amp*s_peaksSD;          % Amplitud máxima = +8 SD
    
    [v_peaks, v_peaksind] = findpeaks (v_ecginteg_temp,... 
 'MinPeakHeight', s_minThreshold,'MinPeakDistance',s_mindursam);% Detección picos 
    v_peaksind = v_peaksind(v_peaks <= s_maxThreshold );        % Depuración, de picos mayores

    v_peaksind_Cell{1,i} = v_peaksind;                          % Incluir vector temporal en la célula
    v_peaks_Cell{1,i}    = v_peaks; 
    
    clear v_ecginteg_temp v_peaks s_peaksMean s_peaksSD ...
        s_minThreshold s_maxThreshold v_peaksind
end                              % Detección de -R-
     % Ajuste por corrimiento producido integración (Verdaderos R)
s_QRSPeakAdjustHalfWinSec = 0.05;                               % Duración QRS (s)
s_QRSPeakAdjustHalfWinSam = ...                                 
    round(s_QRSPeakAdjustHalfWinSec*s_FsHz);                    % Duración QRS (muestras)
for i = 1:length(v_ECGSig_Cell)
  v_peaksind(1,:) = v_peaksind_Cell{1,i};                            % Variables temporales 
  
  v_ecgfilt  = v_ECGSig_Cell{1,i};                              % 
    for s_count = 1:numel(v_peaksind)
        s_ind = v_peaksind(s_count);                            % Variable temporal para guardar 
        s_firstind = s_ind -(s_QRSPeakAdjustHalfWinSam);
        s_lastind  = s_ind +(s_QRSPeakAdjustHalfWinSam);
            if s_firstind < 1
                s_firstind = 1; 
            end
            if s_lastind > max(v_peaksind)
                s_lastind = max(v_peaksind); 
            end         
        [~, s_ind] = findpeaks(v_ecgfilt(s_firstind:s_lastind),...
          'SortStr', 'descend');                            % Eliminar falsos positivos
             if isempty(s_ind)
                continue
             end 
        s_ind = s_ind(1);
        v_peaksind(s_count) = s_firstind + s_ind -1; % Correción desfase 
        clear s_ind s_firstind s_lastind
    end                        % Detección de verdaderos positivos
  v_peaksind_Cell{1,i} = []; 
  v_peaksind_Cell{1,i} = v_peaksind; 
  clear v_peaksind v_ecgfilt
end                              % Corrección de desfase
%% Tacograma - bpm 
v_Taco_cell{1,length(v_ECGSig_Cell)} = [];                      % Preasignación de células
v_bpm_cell{1,length(v_ECGSig_Cell)}  = [];                      % 
for i = 1: length(v_ECGSig_Cell)
   v_peaksind       = v_peaksind_Cell{1,i};
   v_Taco           = diff(v_peaksind)./ s_FsHz;                % Intervalo R-R instantáneo
   bpm              = 60./v_Taco;                               % Latidos por minuto instantaneos 
   v_Taco_cell{1,i} = v_Taco;                               
   v_bpm_cell{1,i}  = bpm; 
   clear v_Taco bpm v_peaksind
end                             % Tacograma y bpm

%% Estadísticos de tacograma Y bpm
% (Normal to normal) = Normalmente se excluyen de la distribución
%  +1SD +8SD. 

% 1.0 Remuestrear tacograma & BPM                               % En este paso se pierde 
 
v_Taco_cell_res {1,length(v_ECGSig_Cell)} = [];                 % Preasignación de célula
v_bpm_cell_res {1,length(v_ECGSig_Cell)}  = [];    
for i = 1:length(v_ECGSig_Cell)
 v_Taco_temp = v_Taco_cell{1,i};
 peak_temp = v_peaksind_Cell{1,i};
 [v_Taco_temp_res,~] = resample(v_Taco_temp,peak_temp(2:end),1,3,1,'linear');
 v_Taco_cell_res{1,i} = v_Taco_temp_res;

 clear v_Taco_temp_res v_Taco_temp peak_temp
 end                             % algo de precisión
 for i = 1:length(v_ECGSig_Cell)
     v_bpm_temp = v_bpm_cell{1,i};
 peak_temp = v_peaksind_Cell{1,i};
 [v_bpm_temp_res,~] = resample(v_bpm_temp,peak_temp(2:end),1,3,1,'linear');
 v_bpm_cell_res{1,i} = v_bpm_temp_res;

 clear v_bpm_temp_res v_bpm_temp peak_temp
 end 
%%% Parámetros para el cálculo de los estadísticos 
s_window = 300;                                                 % Duración ventana (s), 5min
s_windowSam = ...                                 
    round(s_window*s_FsHz);                                     % Duración ventana (muestras)
s_corr = 5;                                                     % Corrimiento (s)
s_corrSam = ...                                 
    round(s_corr*s_FsHz);                                       % Corrimiento(muestras)
% 1.1 Media de los intervalos NN (Normal to normal)
v_meanNN_Cell{1,7} = [];                                        
for i = 1:length(v_ECGSig_Cell)
  v_Taco_temp(1,:) = v_Taco_cell_res{1,i};                      % Variable temporal
  meanNN_temp = zeros (1,numel(v_Taco_temp)/s_corrSam);
  for s_count = 1:numel(v_Taco_temp)/s_corrSam
        s_ind = (s_count-1)*s_corrSam;                          % Factor de corrimiento
        s_firstind = round(s_ind) ;                                    % Punto inicial
            if s_count == 1
                s_firstind = 1; 
            end
        s_lastind  = round(s_ind +(s_windowSam));                      % Punto final 
            if s_lastind > length(v_Taco_temp)
                s_lastind = length(v_Taco_temp); 
            end         
            if s_firstind > length(v_Taco_temp)
                s_firstind = length(v_Taco_temp); 
            end  
     meanNN_temp(1,s_count) = mean(v_Taco_temp(1,s_firstind:s_lastind)); % Media N-N               
  end
  
  v_meanNN_Cell{1,i} = meanNN_temp;
  clear meanNN_temp s_ind s_firstind s_lastind v_Taco_temp
end                              % Media de los intervalos                  
% 1.2 Desviación estándar de intervalos NN (Normal to normal)
v_SDNN_Cell{1,7} = [];  
for i = 1:length(v_ECGSig_Cell)
  v_Taco_temp(1,:) = v_Taco_cell_res{1,i};                      % Variable temporal
  v_SDNN_temp = zeros (1,numel(v_Taco_temp)/s_corrSam);  
  for s_count = 1:numel(v_Taco_temp)/s_corrSam
        s_ind = (s_count-1)*s_corrSam;                          % Factor de corrimiento
        s_firstind = round(s_ind) ;                                    % Punto inicial
            if s_count == 1
                s_firstind = 1; 
            end
        s_lastind  = round(s_ind +(s_windowSam));                      % Punto final 
            if s_lastind > length(v_Taco_temp)
                s_lastind = length(v_Taco_temp); 
            end         
            if s_firstind > length(v_Taco_temp)
                s_firstind = length(v_Taco_temp); 
            end  
     v_SDNN_temp(1,s_count) = std(v_Taco_temp(1,s_firstind:s_lastind)); % SD N-N               
  end
  
  v_SDNN_Cell{1,i} = v_SDNN_temp;
  clear v_SDNN_temp s_ind s_firstind s_lastind v_Taco_temp
end
% NN50 = Número de pares adyacentes de intervalos NN que difieren de más de 50 ms
NN50_cell{1,7} = []; 
for i = 1:length(v_ECGSig_Cell)
  v_Taco_temp = v_Taco_cell{1, i}  ;
  NN50_temp = 0;
  for ind = 1:  length(v_Taco_temp)
    if v_Taco_temp(ind) > 0.5
        NN50_temp = NN50_temp + 1;
    end
  end
  NN50_cell{1,i} = NN50_temp;
  clear NN50_temp v_Taco_temp
end
% pNN50 = NN50 dividido por el número total de intervalos NN
pNN50_cell{1,7} = [];
for i = 1:length(v_ECGSig_Cell)
  NNinttot = length(v_Taco_cell{1, i});
  NN50_temp = (NN50_cell{1,i});
  pNN50_temp = NN50_temp/NNinttot;
  pNN50_cell{1,i} = pNN50_temp; % proporción de NN > 50ms
  clear NNinttot v_Taco_temp pNN50_temp
end
% MEANNNBPM = Media de los intervalos NN en latidos por minutos instantáneos.
v_meanbpmNN_Cell{1,7} = [];       
for i = 1:length(v_ECGSig_Cell)
  v_bpm_temp(1,:) = v_bpm_cell_res{1,i};                      % Variable temporal
  meanbpmNN_temp = zeros (1,numel(v_bpm_temp)/s_corrSam);
  for s_count = 1:numel(v_bpm_temp)/s_corrSam
        s_ind = (s_count-1)*s_corrSam;                          % Factor de corrimiento
        s_firstind = round(s_ind) ;                                    % Punto inicial
            if s_count == 1
                s_firstind = 1; 
            end
        s_lastind  = round(s_ind +(s_windowSam));                      % Punto final 
            if s_lastind > length(v_bpm_temp)
                s_lastind = length(v_bpm_temp); 
            end         
            if s_firstind > length(v_bpm_temp)
                s_firstind = length(v_bpm_temp); 
            end  
     meanbpmNN_temp(1,s_count) = mean(v_bpm_temp(1,s_firstind:s_lastind)); % Media N-N               
  end
  
  v_meanbpmNN_Cell{1,i} = meanbpmNN_temp;
  clear meanbpmNN_temp s_ind s_firstind s_lastind v_bpm_temp
end 
% SDNNBPM = Desviación estándar de intervalos NN en latidos por minutos instantáneos.
v_SDNNbpm_Cell{1,7} = [];                                        
for i = 1:length(v_ECGSig_Cell)
  v_bpm_temp(1,:) = v_bpm_cell_res{1,i};                      % Variable temporal
  v_SDNNbpm_temp = zeros(1,numel(v_bpm_temp)/s_corrSam);
  for s_count = 1:numel(v_bpm_temp)/s_corrSam
        s_ind = (s_count-1)*s_corrSam;                          % Factor de corrimiento
        s_firstind = round(s_ind) ;                                    % Punto inicial
            if s_count == 1
                s_firstind = 1; 
            end
        s_lastind  = round(s_ind +(s_windowSam));                      % Punto final 
            if s_lastind > length(v_bpm_temp)
                s_lastind = length(v_bpm_temp); 
            end         
            if s_firstind > length(v_bpm_temp)
                s_firstind = length(v_bpm_temp); 
            end  
     v_SDNNbpm_temp(1,s_count) = std(v_bpm_temp(1,s_firstind:s_lastind)); % SD N-N               
  end
  
  v_SDNNbpm_Cell{1,i} = v_SDNNbpm_temp;
  clear v_SDNNbpm_temp s_ind s_firstind s_lastind v_bpm_temp
end

%% Estadísticos para Fase del sueño:

for suj = 1:length (v_ECGSig_Cell)
    v_meanNN_temp  = v_meanNN_Cell{1,suj};
    v_HypTime      = (0:30:30*numel(v_HypCode_Cell{1,suj}))';
    v_HypTime(end) = [];
    Hyp_time_score = cat(2,v_HypCode_Cell{1,suj},v_HypTime);
    
for k = 1:length(v_HypCodeLabels)
    x = 1;                                                          % Variables de posición en la matriz
    y = 1; 

for i = 1:length (Hyp_time_score)-1
    if Hyp_time_score (i,1) == k 
        Temporal_times(x,y) = Hyp_time_score(i,2); 
        y = y+1;
    elseif Hyp_time_score (i,1) && Hyp_time_score (i+1,1) ~= k 
    else
        x = x+1;
        y = 1;
    end
end                         
Hyp_seg_times{1, k} = Temporal_times;
clear Temporal_times
end                            % Matriz con momentos de WAKE
clear i k x y v_HypTime

temp_vmeanNN = v_meanNN_Cell{1,sujeto}; 
for k = 1:length(v_HypCodeLabels)
    Hyp_seg_temp = Hyp_seg_times{1,k}; 
    for row = 1: size(Hyp_seg_temp,1)
       for col = 1:size(Hyp_seg_temp,2)
        if Hyp_seg_temp(row,col) == 0 
            Hyp_seg_temp(row,col) = NaN(1);                     % Excluir ceros
        end 
       end 
    end 
    clear row col
    Temp_vmeanNN_sleep_cat = [];
   for i = 1:size(Hyp_seg_temp,1)
        min_lim = min(Hyp_seg_temp(i,:)); 
        max_lim = max(Hyp_seg_temp(i,:)); 
        if isnan(min_lim)
            continue
        end 
        if max_lim*0.2 > length(temp_vmeanNN)
            max_lim = length(temp_vmeanNN)/0.2;
        end 
        Temp_vmeanNN_sleep = temp_vmeanNN (min_lim*0.2:max_lim*0.2); % Segmentar por 
        Temp_vmeanNN_sleep_cat = cat(2,Temp_vmeanNN_sleep_cat,Temp_vmeanNN_sleep);
     clear Temp_vmeanNN_sleep min_lim max_lim
   end 
   MEAN_Taco_by_Hyp(suj,k)= mean(Temp_vmeanNN_sleep_cat); 
   clear Temp_vmeanNN_sleep_cat Hyp_seg_temp
end                            % Media para cada fase

temp_v_SDNN = v_SDNN_Cell{1,sujeto}; 
for k = 1:length(v_HypCodeLabels)
    Hyp_seg_temp = Hyp_seg_times{1,k}; 
    for row = 1: size(Hyp_seg_temp,1)
       for col = 1:size(Hyp_seg_temp,2)
        if Hyp_seg_temp(row,col) == 0 
            Hyp_seg_temp(row,col) = NaN(1);                     % Excluir ceros
        end 
       end 
    end 
    clear row col
    Temp_v_SDNN_sleep_cat = [];
   for i = 1:size(Hyp_seg_temp,1)
        min_lim = min(Hyp_seg_temp(i,:)); 
        max_lim = max(Hyp_seg_temp(i,:)); 
        if isnan(min_lim)
            continue
        end 
        if max_lim*0.2 > length(temp_v_SDNN)
            max_lim = length(temp_v_SDNN)/0.2;
        end 
        Temp_v_SDNN_sleep = temp_v_SDNN (min_lim*0.2:max_lim*0.2); % Segmentar por 
        Temp_v_SDNN_sleep_cat = cat(2,Temp_v_SDNN_sleep_cat,Temp_v_SDNN_sleep);
     clear Temp_v_SDNN_sleep min_lim max_lim
   end 
   std_Taco_by_Hyp(suj,k)= std(Temp_v_SDNN_sleep_cat); 
   clear Temp_v_SDNN_sleep_cat Hyp_seg_temp
end                            % SD para cada fase
end

%% Gráficas de tacograma con hipnograma
% Determinar sujeto %%%%
    sujeto = 4;
%%%%%%%%%%%%%%%%%%%%%%%%
v_time        = (0:numel(v_ECGSig_Cell{1,sujeto})-1)./s_FsHz;   % Variables temporales 
v_PeaksInd    = v_peaksind_Cell{1,sujeto};                      %
v_Taco        = v_Taco_cell{1,sujeto};                          %
v_ECGFilt     = v_ECGSig_Cell{1,sujeto};                        %
v_HypCode     = v_HypCode_Cell{1,sujeto};
v_HypTime     = (0:30:30*numel(v_HypCode))';
v_HypTime(end)= [];


figure
subplot(321)                                                    
plot(v_time, v_ECGFilt/1000)                                    % Plot del ECG 
hold on
%plot((v_PeaksInd/s_FsHz), v_ECGFilt(v_PeaksInd), 'o')           % Picos hallados
ylabel('mV');
title('Electrocardiograma');
xlim([v_time(1) v_time(end)])
ylim([min(v_ECGFilt/1000) max(v_ECGFilt/1000)])
grid
xlabel('Tiempo (s)')


subplot(222)
plot(v_time, v_ECGFilt/1000) 
hold on
plot((v_PeaksInd/s_FsHz), v_ECGFilt(v_PeaksInd)/1000, 'o')           % Picos hallados
ylabel('mV');
xlim([32.5 33.5])
ylim([min(v_ECGFilt)/1000 max(v_ECGFilt)/1000])
grid
xlabel('Tiempo (s)')

subplot(323)
hold off

plot(v_time(v_PeaksInd(2:end)), v_Taco)                         % Plot del tacograma (general)
grid
ylabel('R-R (sec.)');
title('Tacograma');
xlim([v_time(1) v_time(end)])
ylim([min(v_Taco) max(v_Taco)])
%media
hold on
v_time_2 = (0:numel(v_meanNN_Cell{1,sujeto})-1)./0.2;           % vector de tiempo (Fs=0.2Hz)
plot(v_time_2,v_meanNN_Cell{1,sujeto})                          % plot media   
%std
hold on
v_time_2 = (0:numel(v_meanNN_Cell{1,sujeto})-1)./0.2;           % vector de tiempo (Fs=0.2Hz)
plot(v_time_2,v_meanNN_Cell{1,sujeto}+ v_SDNN_Cell{1,sujeto}) % plot std
hold on
plot(v_time_2,v_meanNN_Cell{1,sujeto}- v_SDNN_Cell{1,sujeto}, 'color', [0.9290 0.6940 0.1250]) % plot std
xlabel('Tiempo (s)')




subplot(426)
plot(v_time, v_ECGFilt/1000) 
hold on
plot((v_PeaksInd/s_FsHz), v_ECGFilt(v_PeaksInd)/1000, 'o')           % Picos hallados
ylabel('mV');
xlim([30 40])
ylim([min(v_ECGFilt)/1000 max(v_ECGFilt)/1000])
grid

subplot(428)
plot(v_time(v_PeaksInd(2:end)), v_Taco)                         % Plot del tacograma (detalle)
grid
ylabel('R-R (sec.)');
title('Tacograma');
xlim([30 40])
ylim([min(v_Taco) max(v_Taco)])

xlabel('Tiempo (s)')

v_HypTime            = v_HypTime + 15;                          % Parametros hipgnograma
v_Ind                = cell2mat(v_HypCodeLabels(:, 3));
v_HypCodeTicksLabels = cell(1, numel(v_Ind));
v_HypCodeTicks       = zeros(1, numel(v_Ind));
for s_Count = 1:numel(v_Ind)
    v_HypCodeTicks(s_Count)       = str2double(v_Ind(s_Count));
    v_HypCodeTicksLabels{s_Count} = v_HypCodeLabels{s_Count, 2};
end
[~, v_Ind]           = sort(v_HypCodeTicks);
v_HypCodeTicks       = v_HypCodeTicks(v_Ind);
v_HypCodeTicksLabels = v_HypCodeTicksLabels(v_Ind);
       subplot(325)               
plot(v_HypTime, v_HypCode)                                      % Plot del hipnograma
xlim([v_time(1) v_time(end)])
set(subplot(325), 'ytick', v_HypCodeTicks, ...
    'yticklabel', v_HypCodeTicksLabels);
ylabel('Estados de sueño');
title('Hipnograma');
grid
xlabel('Tiempo (s)')

clear v_time  v_PeaksInd   v_ECGFilt  v_HypCode  v_HypTime ...
    v_Axes sujeto v_Ind v_HypCodeTicksLabels v_HypCodeTicks

