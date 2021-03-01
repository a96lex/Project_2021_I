set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energies.png"
set title "Energies respecte temps"
set xlabel "Temps"
set ylabel "Energia"
set offsets 0, 0, 50, 0
plot "./results/thermodynamics.dat" u 1:2 w l t "Energia cinètica","" u 1:3 w l t "Energia potencial","" u 1:4 w l t "Energia total"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Pressio.png"
set title "Pressió respecte temps"
set xlabel "Temps"
set ylabel "Pressió"
set offsets 0, 0, 0, 0
plot "./results/mean_press.dat" u 1:2 w l t "Pressió"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Temperatura.png"
set title "Temperatura respecte temps"
set xlabel "Temps"
set ylabel "Temperatura"
set offsets 0, 0, 0, 0
plot "./results/thermodynamics.dat" u 1:5 w l t "Temperatura"

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energia_cinetica_bins.png"
set title "Variància de la energia cinètica respecte número de bins"
set xlabel "Número de bins"
set ylabel "Variància de la energia cinètica"
set offsets 1, 1, 0, 0
unset key
plot "./results/ekinBIN.dat" u 1:3 w p ls 7 lc rgb "blue" ps 2

set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energia_potencial_bins.png"
set title "Variància de la energia potencial respecte número de bins"
set xlabel "Número de bins"
set ylabel "Variància de la energia potencial"
set offsets 1, 1, 0.1, 0
unset key
plot "./results/epotBIN.dat" u 1:3 w p ls 7 lc rgb "blue" ps 2 