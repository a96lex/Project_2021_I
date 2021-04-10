# Author: Alex

# Configuració global
data_directory = "./results/"
data_directory_dim = "./results/dimensionalized/"
plots_directory = "./results/plots/"
plots_directory_dim = plots_directory."dimensionalized/"
set fit quiet
set fit logfile '/dev/null'
set terminal png size 800,600


set output plots_directory."Energies.png"
set title "Energies respecte temps"
set xlabel "Temps"
set ylabel "Energia"
set offsets 0, 0, 50, 0
set key bmargin center horizontal spacing 3
plot data_directory."/thermodynamics.dat" u 1:2 w l t "Energia cinètica","" u 1:3 w l t "Energia potencial","" u 1:4 w l t "Energia total"


set output plots_directory."Pressio.png"
set title "Pressió respecte temps"
set xlabel "Temps"
set ylabel "Pressió"
set offsets 0, 0, 0, 0
unset key
plot data_directory."thermodynamics.dat" u 1:6 w l 


set output plots_directory."Temperatura.png"
set title "Temperatura respecte temps"
set xlabel "Temps"
set ylabel "Temperatura"
set offsets 0, 0, 0, 0
unset key
plot data_directory."thermodynamics.dat" u 1:5 w l 


set output plots_directory."Energia_cinetica_bins.png"
set title "Binning de la energia cinètica"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia cinètica"
set offsets 1, 1, 0, 0
unset key
plot data_directory."ekinBIN.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2


set output plots_directory."Energia_potencial_bins.png"
set title "Binning de la energia potencial"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia potencial"
set offsets 1, 1, 0.1, 0
unset key
plot data_directory."epotBIN.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2 


set output plots_directory."Distribucio_radial.png"
set title "Distribució radial del sistema"
set xlabel "Distància reduida"
set ylabel "Frequència"
set yrange[0:]
unset key
unset offsets
plot data_directory."radial_distribution.dat" u 1:2:3 w yerrorbars lc rgb "light-blue", \
'' u 1:2 w l lc rgb "blue";


set output plots_directory."Coef_dif_x.png"
set title "Coeficient de difusió X"
set xlabel "2t"
set ylabel "Variància X"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory."diffcoeff.dat" u 1:2 via a,b
plot data_directory."diffcoeff.dat" u 1:2 w l, lin(x) w l lc rgb "blue"


set output plots_directory."Coef_dif_y.png"
set title "Coeficient de difusió Y"
set xlabel "2t"
set ylabel "Variància Y"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory."diffcoeff.dat" u 1:3 via a,b
plot data_directory."diffcoeff.dat" u 1:3 w l, lin(x) w l lc rgb "blue"


set output plots_directory."Coef_dif_z.png"
set title "Coeficient de difusió Z"
set xlabel "2t"
set ylabel "Variància Z"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory."diffcoeff.dat" u 1:4 via a,b
plot data_directory."diffcoeff.dat" u 1:4 w l, lin(x) w l lc rgb "blue"


set output plots_directory."Correlation_energy.png"
set title "Funció Autocorrelació Energia Total"
set xlabel "Time Lag"
set ylabel "Autocorrelació"
unset key
plot data_directory."correlation_energy.dat" w l



# dimensionalized plots

set output plots_directory_dim."Energies.png"
set title "Energies respecte temps"
set xlabel "Temps [s]"
set ylabel "Energia [J/mol]"
set offsets 0, 0, 50, 0
set key bmargin center horizontal spacing 3
plot data_directory_dim."thermodynamics_dim.dat" u 1:2 w l t "Energia cinètica","" u 1:3 w l t "Energia potencial","" u 1:4 w l t "Energia total"


set output plots_directory_dim."Pressio.png"
set title "Pressió respecte temps"
set xlabel "Temps [s]"
set ylabel "Pressió [pa]"
set offsets 0, 0, 0, 0
unset key
plot data_directory_dim."thermodynamics_dim.dat" u 1:6 w l 


set output plots_directory_dim."Temperatura.png"
set title "Temperatura respecte temps"
set xlabel "Temps [s]"
set ylabel "Temperatura [K]"
set offsets 0, 0, 0, 0
unset key
plot data_directory_dim."thermodynamics_dim.dat" u 1:5 w l 


set output plots_directory_dim."Energia_cinetica_bins.png"
set title "Binning de la energia cinètica"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia cinètica [J/mol]"
set offsets 1, 1, 0, 0
unset key
plot data_directory_dim."ekinBIN_dim.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2


set output plots_directory_dim."Energia_potencial_bins.png"
set title "Binning de la energia potencial"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia potencial [J/mol]"
set offsets 1, 1, 0.1, 0
unset key
plot data_directory_dim."epotBIN_dim.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2 


set output plots_directory_dim."Distribucio_radial.png"
set title "Distribució radial del sistema"
set xlabel "Distància [{\305}]"
set ylabel "Frequència"
set yrange[0:]
unset key
unset offsets
plot data_directory_dim."radial_distribution_dim.dat" u 1:2:3 w yerrorbars lc rgb "light-blue", \
'' u 1:2 w l lc rgb "blue";


set output plots_directory_dim."Coef_dif_x.png"
set title "Coeficient de difusió X"
set xlabel "Temps [s]"
set ylabel "Variància X [{\305}]"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory_dim."diffcoeff_dim.dat" u 1:2 via a,b
plot data_directory_dim."diffcoeff_dim.dat" u 1:2 w l, lin(x) w l lc rgb "blue"


set output plots_directory_dim."Coef_dif_y.png"
set title "Coeficient de difusió Y"
set xlabel "Temps [s]"
set ylabel "Variància Y [{\305}]"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory_dim."diffcoeff_dim.dat" u 1:3 via a,b
plot data_directory_dim."diffcoeff_dim.dat" u 1:3 w l, lin(x) w l lc rgb "blue"


set output plots_directory_dim."Coef_dif_z.png"
set title "Coeficient de difusió Z"
set xlabel "Temps [s]"
set ylabel "Variància Z [{\305}]"
set yrange[*:*]
unset key
lin(x)=a*x+b
fit lin(x) data_directory_dim."diffcoeff_dim.dat" u 1:4 via a,b
plot data_directory_dim."diffcoeff_dim.dat" u 1:4 w l, lin(x) w l lc rgb "blue"


set output plots_directory_dim."Correlation_energy.png"
set title "Funció Autocorrelació Energia Total"
set xlabel "Time Lag"
set ylabel "Autocorrelació"
unset key
plot data_directory_dim."correlation_energy_dim.dat" w l