# Author: Alex

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energies.png"
set title "Energies respecte temps"
set xlabel "Temps"
set ylabel "Energia"
set offsets 0, 0, 50, 0
set key bmargin center horizontal spacing 3
plot "./results/thermodynamics.dat" u 1:2 w l t "Energia cinètica","" u 1:3 w l t "Energia potencial","" u 1:4 w l t "Energia total"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Pressio.png"
set title "Pressió respecte temps"
set xlabel "Temps"
set ylabel "Pressió"
set offsets 0, 0, 0, 0
unset key
plot "./results/thermodynamics.dat" u 1:7 w l 

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Temperatura.png"
set title "Temperatura respecte temps"
set xlabel "Temps"
set ylabel "Temperatura"
set offsets 0, 0, 0, 0
unset key
plot "./results/thermodynamics.dat" u 1:5 w l 

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energia_cinetica_bins.png"
set title "Binning de la energia cinètica"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia cinètica"
set offsets 1, 1, 0, 0
unset key
plot "./results/ekinBIN.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energia_potencial_bins.png"
set title "Binning de la energia potencial"
set xlabel "Número de bins"
set ylabel "Desviació quadràtica de la energia potencial"
set offsets 1, 1, 0.1, 0
unset key
plot "./results/epotBIN.dat" u 1:3 w lp ls 7 lc rgb "blue" ps 2 

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Distribucio_radial.png"
set title "Distribució radial del sistema"
set xlabel "Distància reduida"
set ylabel "Frequència"
set yrange[0:]
unset key
unset offsets
plot "./results/radial_distribution.dat" u 1:3:4 w yerrorbars lc rgb "light-blue", \
'' u 1:3 w l lc rgb "blue";

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Coef_dif_x.png"
set title "Coeficient de difusió X"
set xlabel "2t"
set ylabel "Variància X"
unset key
lin(x)=a*x+b
fit lin(x) "./results/diffcoeff.dat" u 1:2 via a,b
plot "./results/diffcoeff.dat" u 1:2 w l, lin(x) w l lc rgb "blue"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Coef_dif_y.png"
set title "Coeficient de difusió Y"
set xlabel "2t"
set ylabel "Variància Y"
unset key
lin(x)=a*x+b
fit lin(x) "./results/diffcoeff.dat" u 1:3 via a,b
plot "./results/diffcoeff.dat" u 1:3 w l, lin(x) w l lc rgb "blue"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Coef_dif_z.png"
set title "Coeficient de difusió Z"
set xlabel "2t"
set ylabel "Variància Z"
unset key
lin(x)=a*x+b
fit lin(x) "./results/diffcoeff.dat" u 1:4 via a,b
plot "./results/diffcoeff.dat" u 1:4 w l, lin(x) w l lc rgb "blue"

