
function optionPrice = EuropeanOptionAmericanBarrier(F0, K, KO, B, TTM, sigma, d)

    % Funzione per calcolare il prezzo di una Call Up-and-Out (Barriera Americana)
    % utilizzando il Principio di Riflessione come mostrato nell'assignment.
    
    % 1. Ricaviamo il tasso r dal fattore di sconto B
    % Assumendo che B = exp(-r * TTM)
    r = -log(B) / TTM;
    S_curr = F0.*B.*exp(d*TTM);
    
    % 2. Controllo barriera (Knock-Out immediato)
    if S_curr >= KO
        optionPrice = 0;
    else
        % 3. Calcolo del Forward price corrente
        F0_curr = S_curr * exp(-d * TTM) / B;
        
        % 4. Prezzo dell'opzione Europea standard
        % NOTA: Assicurati che la funzione 'OptionPrice_closed_eur_barrier' 
        % sia presente nella tua stessa cartella di lavoro
        P_u = EuropeanKOCall_ClosedFormula(F0_curr, K, KO, B, TTM, sigma);
        
        % 5. Calcolo del Forward price "riflesso" (Mirrored)
        F1_curr = (KO^2 / S_curr) * exp(-d * TTM) / B;
        
        % 6. Prezzo dell'opzione Europea "riflessa"
        P_u1 = EuropeanKOCall_ClosedFormula(F1_curr, K, KO, B, TTM, sigma);
        
        % 7. Principio di Riflessione per il calcolo del prezzo finale
        esponente = 2 * (r - d - 0.5 * sigma^2) / sigma^2;
        optionPrice = P_u - (KO / S_curr)^esponente * P_u1;
    end
end