%% REPRESENTACIÓN DE PROPIEDADES DE SOFA

% Leemos los valores del archivo txt generado de Sofa y lo metemos en una
% tabla
T = readtable('properties_0-28-56_x.txt');

% Desglosamos los valores en sus respectivas componentes
[vector_tiempo,posicion_nodo_0_dcha,posicion_nodo_28_central,posicion_nodo_56_izq] = readvars('properties_0-28-56_x.txt');

 % Expresion que ayuda a dividir la celda cada que encontramos un espacio en el vector
expression = ' ';

% Dividimos la cadena en la misma variable pero con celdas separadas según
% hemos indicado en la variable expresión. Es decir que ashora por cada
% fila tendremos a su vez una celda dividida en 7 celdas más pequeñas
posicion_nodo_0_dcha_split = regexp(posicion_nodo_0_dcha,expression,'split');
posicion_nodo_28_central_split = regexp(posicion_nodo_28_central,expression,'split');
posicion_nodo_56_izq_split = regexp(posicion_nodo_56_izq,expression,'split');

for i=1:1032
    posicion_nodo_0_dcha_split_ = posicion_nodo_0_dcha_split{i,1}; % Variable that has the numerical value of latitude.
    % posiciones(i,7) = posicion_nodo_0_dcha_split_
    for j=1:7
        posicion_nodo_0_dcha_split_num = cell2mat(posicion_nodo_0_dcha_split_(1,j)); % Convert from cell to numerical value.
    end
end