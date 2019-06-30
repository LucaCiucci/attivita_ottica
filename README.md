# attivita_ottica

## note
I sorgenti sono solo per dare un'idea degli algoritmi. Non possono essere compilati così come sono, forse funzionano includendo gli header nella git LC e SFML. LC è diversa da quella utilizzata nel vero progetto, potrebbe non essere compatibile ma con delle modifiche dovrebbe andare. Comunque l'eseguibile funzionante per Windows 64 è in /bin.  
Mi scuso per il disastro in questo codice, ma a forza di modifiche e cose provvisorie che ho fatto durante tutta l'esperienza è diventato un disastro, specie il main che ho usato anche per fare prove di cose tipo F-test. comunque le cose interessanti stanno in haloFitter.inl.  
## utilizzo
lanciare l'eseguibile. Dopo un po' di linee inutili chiede il nome dell'immagine da utilizzare, scrivendo 1 o 2/3/4/5/6 prende una delle immagini fornite.  
0..2 per il canale da utilizzare (R/G/B).  
iIl parametro successivo lasciare a 0.  
Poi il numero di punti in cui suddividere l'angolo giro.  
Si apre una finestra. Con cliccare su alcuni punti per una circonferenza approssimativa. inoltre i seguenti comandi:  
f2 per attivare/disattivare la circonferenza  
hold R/G/B per visualizzare solo le componenti desiderate  
+/- per aumentare la luminosità  
rotellina +/- per zoom  
rotellina click per pan  
  
premendo ESC poi parte la fase successiva. L'utente deve premere Y/N per accettare o meno il punto scelto.
