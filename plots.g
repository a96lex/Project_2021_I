set terminal png size 800,600
system "mkdir -p ./results/plots"
set output "./results/plots/Energies.png"
set title "Plot test"
plot "./results/thermodynamics.dat" u 1:2 w l t "Energia cinètica","" u 1:3 w l t "Energia potencial","" u 1:4 w l t "Energia total"

