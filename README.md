Información sobre los códigos:

# Atenuación_Ha.py

Este código toma datos de un archivo del Millenium para la línea H_alpha y dibuja, para un redshift z dado, la función de luminosidad sin atenuar, la atenuada por Millenium, la atenuada por Calzetti y la de por Cardelli. Se asume un caja de lado 200 Mpc/h. Necesita usar Att_Calzetti y Att_Cardelli.
Resultado: figura 3 del TFG


# Atenuación_Hb.py

Este código toma datos de un archivo del Millenium para la línea H_beta y dibuja, para un redshift z dado, la función de luminosidad sin atenuar, la atenuada por Millenium, la atenuada por Calzetti y la de por Cardelli y los datos experimentales de Comparat 2016. Se asume un caja de lado 200 Mpc/h. Necesita usar Att_Calzetti y Att_Cardelli.
Resultado: figura 4 del TFG


# Atenuación_OII.py

Este código toma datos de un archivo del Millenium para el doblete OII y dibuja, para un redshift z dado, la función de luminosidad sin atenuar, la atenuada por Millenium, la atenuada por Calzetti y la de por Cardelli y los datos experimentales de Comparat 2016. Se asume un caja de lado 200 Mpc/h. Necesita usar Att_Calzetti y Att_Cardelli.
Resultado: figura 5 del TFG


# Atenuación_OII_z.py

Este código toma datos de un archivo del Millenium para el doblete OII y calcula, para un redshift z dado, la función de luminosidad sin atenuar, la atenuada por Millenium, la atenuada por Calzetti y la de por Cardelli. Se crean 3 gráficas según el tipo de atenuación: Cardelli, Calzetti y Millenium. Se asume un caja de lado 200 Mpc/h. Necesita usar Att_Calzetti y Att_Cardelli.
Resultado: figuras 6,7 y 8 del TFG


# AttOII_Saito.py

Este código toma datos de un archivo del Millenium para el doblete OII y dibuja, para un redshift z dado, la función de luminosidad sin atenuar, la atenuada por Calzetti, la de por Calzetti con el factor de Saito y los datos experimentales de Comparat 2016. Se asume un caja de lado 200 Mpc/h. Necesita usar Att_Calzetti.
Resultado: figuras 9,10,11 y 12 del TFG


# Att_Calzetti.py

Función de Calzetti para meter atenuación al logaritmo de la luminosidad, se basa en el cálculo de Favole para la profundidad óptica. Necesita usar Cosmology.py y los datos de las galaxias de Millenium

# Att_Cardelli.py

Función de Cardelli para meter atenuación al logaritmo de la luminosidad, se basa en el cálculo de Favole para la profundidad óptica. Necesita usar Cosmology.py y los datos de las galaxias de Millenium

# Cosmology.py

Función auxiliar para definir una Cosmología. Creada por Violeta González Pérez

# plot9_mejorado.py

Función de luminosidad de los modelos de MultiDark con atenuación de Cardelli
Resultado: figura 2 del TFG
