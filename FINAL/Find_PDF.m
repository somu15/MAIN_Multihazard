function [Fr] = Find_PDF(tr,E,PHI)

Ed = Differentiation(tr(2)-tr(1),E);
PHId = Differentiation(tr(2)-tr(1),PHI);

Fr = Ed.*PHI+E.*PHId;

C = trapz(tr,Fr);

Fr = Fr/C;

end