# Simulation of Fluid Droplet using Kawasaki Dynamics

The simulation uses the following hamiltonian:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}&space;=&space;J&space;\sum_{i,j:nn}&space;\sigma_i\sigma_j&space;-&space;K&space;\sum_{i,j:nnn}&space;-&space;\sum_j&space;h_j&space;\sum_{\text{line&space;}j}&space;\sigma_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}&space;=&space;J&space;\sum_{i,j:nn}&space;\sigma_i\sigma_j&space;-&space;K&space;\sum_{i,j:nnn}&space;-&space;\sum_j&space;h_j&space;\sum_{\text{line&space;}j}&space;\sigma_i" title="\mathcal{H} = J \sum_{i,j:nn} \sigma_i\sigma_j - K \sum_{i,j:nnn} - \sum_j h_j \sum_{\text{line }j} \sigma_i" /></a>

It was carried out using Metropolis Monte Carlo with the acceptance probability

<a href="https://www.codecogs.com/eqnedit.php?latex=A(X&space;\to&space;Y)&space;=&space;\min\{1,&space;e^{\beta&space;\Delta&space;E}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A(X&space;\to&space;Y)&space;=&space;\min\{1,&space;e^{\beta&space;\Delta&space;E}\}" title="A(X \to Y) = \min\{1, e^{\beta \Delta E}\}" /></a>
