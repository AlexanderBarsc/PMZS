function [Wo, n] = VygenerujObdelnikovouOkenniFunkci(pocetVzorku, pocatek, konecPozitivnichHodnot)
%UNTITLED Tato funkce vygeneruje obdelnikovou okenni funkci
%   Detailed explanation goes here

konec = pocetVzorku + pocatek;

n = pocatek:1:konec - 1;
Wo = zeros(1, pocetVzorku)

for i = 1:length(n)
    if i > konecPozitivnichHodnot + 1
        Wo(1,i) = 0;
    else
        Wo(1,i) = 1;
    end
end


