Il programma eseguirà i seguenti calcoli:

- Scrivere una funzione che calcoli lo spettro di potenza di una serie numerica x(n) usando il metodo di Welch (WOSA). Tutti e 3 i tipi di finestra visti a lezione (Bartlett, Welch e Hann) devono essere implementati lasciando la scelta di quale utilizzare all’utente.

- Lo spettro di potenza deve essere normalizzato in modo che l’area sotto lo spettro di potenza  sia uguale alla varianza di x(n), inoltre prima di calcolare lo spettro di potenza x(n) deve essere corretto affinché abbia una media nulla e eliminando l’eventuale trend lineare.

- Sperimentare la funzione ottenuta su un insieme dei dati visti a lezione (per es. Ox677_1Ma) e commentare le differenze ottenute usando finestre di diverso tipo e un diverso numero di segmenti.

Il prototipo della funzione sarà simile (ma non necessariamente uguale) al seguente:

[spectrum, freq] = mcwosa(x, nseg, overlap, window) 

x: serie temporale con i dati da analizzare

nseg: numero di segmenti utilizzato per calcolare lo spero di potenza

overlap: sovrapposizione dei segmenti (per es. 0.5 indica una sovrapposizione del 50%).

window: finestra usata per ridurre il problema del leakage.


La funzione deve restituire un array con lo spettro di potenza  (spectrum) e un altro array con la frequenza (freq). Per semplificare possiamo assumere che l’intervallo di campionamento sia unitario.
