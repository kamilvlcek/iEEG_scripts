function [t]=tabs2t(tabs)
% vytvori pole t z pole tabs
% tabs jsou Serial Date Number, t je cas v sekundach od zacatku
t = (tabs-tabs(1)).*24.*3600;
end