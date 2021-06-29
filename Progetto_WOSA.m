% Studente Attarantato Kevin %%%%%%%%%%%%%%%%%%%%%%%%%%
% matricola N. 282785 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corso di modelizazzione geologica %%%%%%%%%%%%%%%%%%%
% Anno accademico 2018-2019 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Progetto sul calcolo del power spectrum e frequenza %

      function[spettro, frequency] = Progetto_WOSA(dati, nseg, overlap, finestra)    % inizio funzione chiamata "Progetto WOSA"
      pkg load signal;                                % caricamento del package "signal"
      varianza = var(dati(:,2),1,1);                  % calcolo la varianza
      p = dati(:,2)- mean(dati(:,2));                 % prendo tutta la seconda colonna - min valore della media seconda colonna      
      p = detrend(p, 1);                              % detrendizzo
      lunp = length(p);                               % carico in lunp(p perchè riferita ai dati) la lunghezza dei dati 
      lunghezza_segmento = lunp/nseg;                 % decido la lunghezza dei segmenti (sulla base del vettore) 
      modulo = mod(lunp, nseg);                       % resto divisione
      sovrapposizione = lunghezza_segmento/overlap;   % la sovrapposizione è uguale alla lungh. del segmento / i segmenti scelti
      while(modulo!=0 || mod(sovrapposizione, 1)!=0)  % ciclo while per evitare di avere perdita di dati agli estremi 
        p(end +1) = 0;                                % aggiungo uno 0 in posizione finale 
        lunp += 1;                                    % aggiorno la lunghezza
        modulo = mod(lunp, nseg);                     % resto divisione per uscire dal while 
        lunghezza_segmento = lunp/nseg;               % la lunghezza del segmento sarà uguale alla lunghezza di p / il numero dei segmenti 
        sovrapposizione = lunghezza_segmento/overlap; % la sovrapposizione sarà uguale alla lunghezza del segmento / la sovrapposizione scelta dall'utente
      endwhile
      
      % separo i segmenti 
      matrice_p_s = zeros(nseg, lunghezza_segmento);        % nseg righe e lunghezza_segmento colonne, matrice_p_s sta per matrice segmenti
      for indice  = 1:nseg;                                 % ciclo for che parte dall'inizio fino al n. di segmenti
        if(i == 0)                                          % per i = o la matrice mi va dall'inizio alla fine del segmento
          matrice_p_s(indice ,:) = p(1:lunghezza_segmento); % prendo gli elementi alla lunghezza giusta (dipende dalla lunghezza del segmento stesso e dalla sovrapposizione) e lo infilo dentro una matrice che è divisa su piu righe e colonne (dove ogni riga e ogni colonna contiene un pezzo segmentato)
        else                                                % altrimenti la matrice andrà dall'indice - 1* la sovrapposizione+1 alla lunghezza del segmento + la sovrapposizione* l'indice -1
          matrice_p_s(indice ,:) = p((indice-1)*sovrapposizione+1:lunghezza_segmento+sovrapposizione*(indice-1));
        endif                                               % fine del if
      endfor                                                % fine del for
	  
    % normalizzo
	  for indice = 1:nseg;                                  % ciclo che va da 1 al numero dei segmenti
		matrice_p_s(indice,:) /= std(matrice_p_s(indice,:));  % divido la matrice creata precedentemente per i segmenti per la deviazione standard della matrice
	  endfor                                                % fine ciclo for

      switch (finestra)                                   % inizio switch per la scelta delle finestre
        case 1                                            % caso 1
          windows = welchwin(length(matrice_p_s(1, :)));  % scegliamo di applicare la finestra welch calcolata su tutta la matrice dei segmenti
        case 2                                            % caso 2
          windows = bartlett(length(matrice_p_s(1, :)));  % scegliamo di applicare la finestra bartlett calcolata su tutta la matrice dei segmenti
        otherwise                                         % caso 3
          windows = hann(length(matrice_p_s(1, :)));      % scegliamo di applicare la finestra hann calcolata su tutta la matrice dei segmenti
        endswitch    
        
      % applicco la finestra scelta
      for indice = 1:nseg                                             % inizio del ciclo for che va da 1 ai numeri di segmenti
        matrice_p_f(indice, :) = matrice_p_s(indice ,:) .* windows';  % creiamo cosi la matrice con i dati "finestrati"
      endfor                                                          % fine del ciclo for


      % calcolo lo spettro dei segmenti  
      for indice = 1:nseg                                                                   % ciclo for che parte da 1 fino al numero di segmenti     
        matrice_p_sp(indice, :) = (abs(fft(matrice_p_f(indice ,:))/lunghezza_segmento)).^2; % creiamo la matrice con i dati "spettrati"
      endfor                                                                                % fine del ciclo for
      
      % calcolo la media degli spettri 
      spettro = zeros(length(matrice_p_sp(1,:)));             % creo una matrice di 0 con lunghezza della matrice "spettrata"
      for indice = 1:nseg                                     % ciclo for che va da 1 al numero di segmenti
        spettro .+= matrice_p_sp(indice, :);                  % lo spettro e uguale alla somma di tutti gli elementi della matrice "spettrata"
      endfor                                                  % fine del ciclo for
      spettro = (spettro ./ nseg);                            % abbiamo che lo spettro è uguale allo spettro precedentemente ottenuto / il numero di segmenti 
      
      spettro = resize(spettro', length(spettro)/2+1, 1);     % uso la resize per "ritagliare" il vettore, ovvero diminuiamo la grandezza del vettore 
      
      %aggiorno la lunghezza del segmento
      lunghezza_segmento = length(spettro);                   % la lunghezza del segmento sarà uguale alla lunghezza dello spettro
      
      %prendo le frequenze fino a 1/delta*2 (frq.nyquist) ed essendo delta unitario sappiamo che arriva fino a 0
      frequency = [0:length(spettro)-1]/length(spettro);      % calcolo la frequenza 
      lung_fr = length(frequency);                            % setto la lunghezza della frequenza
      lung_fr = lung_fr/2 + 1;                                % e la aggiorno dividendola per 2 + 1
      frequency = resize(frequency, 1, lung_fr);              % "ritaglio" con la funzione resize la parte della frequenza che mi interessa

  
endfunction                                                   % fine della funzione