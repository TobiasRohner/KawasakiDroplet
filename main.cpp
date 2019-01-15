#include <cstddef>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <tuple>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>



std::tuple<std::unique_ptr<char[]>, std::vector<long>, std::vector<long>, double> generate_grid(const long width, const long height, const double J, const double K, const double g, const long dropwidth, const long dropheight) {
    auto grid = std::make_unique<char[]>((width+4)*(height+4));
    // Fill the grid
    for (long y = 0 ; y < 2 ; ++y)
	for (long x = 0 ; x < width+4 ; ++x)
	    grid[(width+4)*y + x] = -1;
    for (long y = 2 ; y < dropheight+2 ; ++y)
	for (long x = 0 ; x < width+2-dropwidth ; ++x)
	    grid[(width+4)*y + x] = -1;
    for (long y = 2 ; y < dropheight+2 ; ++y)
	for (long x = width+2-dropwidth ; x < width+4 ; ++x)
	    grid[(width+4)*y + x] = 1;
    for (long y = dropheight+2 ; y < height+4 ; ++y)
	for (long x = 0 ; x < width+4 ; ++x)
	    grid[(width+4)*y + x] = -1;
    // Fill the particles vector
    std::vector<long> particles;
    particles.reserve((dropwidth+2)*dropheight);
    for (long y = 2 ; y < dropheight+2 ; ++y)
	for (long x = width+2-dropwidth ; x < width+2 ; ++x)
	    particles.push_back((width+4)*y + x);
    // Fill the border vector
    std::vector<long> border;
    border.reserve(dropwidth+dropheight);
    for (long y = 2 ; y < dropheight+2 ; ++y)
	border.push_back((width+4)*y + width+1-dropwidth);
    for (long x = width+2-dropwidth ; x < width+2 ; ++x)
	border.push_back((width+4)*(dropheight+2) + x);
    // Compute the energy
    long Enn = 0;
    for (long y = 2 ; y < height+2 ; ++y)
	for (long x = 2 ; x < width+2 ; ++x)
	    Enn += grid[(width+4)*y + x] * (grid[(width+4)*y + x-1] + grid[(width+4)*y + x+1] + grid[(width+4)*(y-1) + x] + grid[(width+4)*(y+1) + x]);
    long Ennn = 0;
    for (long y = 2 ; y < height+2 ; ++y)
	for (long x = 2 ; x < width+2 ; ++x)
	    Ennn += grid[(width+4)*y + x] * (grid[(width+4)*y + x-2] + grid[(width+4)*y + x+2] + grid[(width+4)*(y-1) + x-1] + grid[(width+4)*(y+1) + x+1] + grid[(width+4)*(y-1) + x+1] + grid[(width+4)*(y+1) + x-1] + grid[(width+4)*(y-2) + x] + grid[(width+4)*(y+2) + x]);
    double Eh = 0;
    for (long y = 2 ; y < height+2 ; ++y) {
	double hj = (2+height-y) * g/height;
	long Ehl = 0;
	for (long x = 2 ; x < width+2 ; ++x)
	    Ehl += grid[(width+4)*y + x];
	Eh += hj * Ehl;
    }
    double E = -J*Enn - K*Ennn - Eh;
    return std::make_tuple(std::move(grid), particles, border, E);
}



template<typename RNG>
static double step(const long width, const long height, const double J, const double K, const double g, const double T,
		   char *grid, std::vector<long>& particles, std::vector<long>& border,
		   RNG& rng) {
    for (;;) {
	// Choose two random places to potentially swap
	size_t ppos_idx, bpos_idx;
	long ppos, bpos;
	std::uniform_int_distribution<size_t> dist_p(0, particles.size()-1);
	for (;;) {
	    std::uniform_int_distribution<size_t> dist_b(0, border.size()-1);
	    ppos_idx = dist_p(rng);
	    bpos_idx = dist_b(rng);
	    ppos = particles[ppos_idx];
	    bpos = border[bpos_idx];
	    // Remove the point from the border if it it no border
	    if (grid[bpos-1] != 1 && grid[bpos+1] != 1 && grid[bpos-width-4] != 1 && grid[bpos+width+4] != 1) {
		std::swap(border[bpos_idx], border.back());
		border.pop_back();
	    }
	    else {
		break;
	    }
	}
	// Compute the difference in energy
	const double dEnn = 2*J * (grid[ppos-1] + grid[ppos+1] + grid[ppos-width-4] + grid[ppos+width+4] -
				   grid[bpos-1] - grid[bpos+1] - grid[bpos-width-4] - grid[bpos+width+4]);
	const double dEnnn = 2*K * (grid[ppos-2] + grid[ppos+2] + grid[ppos-width-5] + grid[ppos+width+5] + grid[ppos-width-3] + grid[ppos+width+3] + grid[ppos-2*(width+4)] + grid[ppos+2*(width+4)] -
				    grid[bpos-2] - grid[bpos+2] - grid[bpos-width-5] - grid[bpos+width+5] - grid[bpos-width-3] - grid[bpos+width+3] - grid[bpos-2*(width+4)] - grid[bpos+2*(width+4)]);
	const long hp = 2+height - ppos/(width+4);
	const long hb = 2+height - bpos/(width+4);
	const double dEh = 2*g/height * (hb - hp);
	const double dE = dEnn + dEnnn + dEh;
	// Accept or reject the swap
	std::uniform_real_distribution<double> dist_acc(0, 1);
	if (dist_acc(rng) <= std::min(1., std::exp(-dE/T))) {
	    // Swap the particle in the grid
	    std::swap(grid[ppos], grid[bpos]);
	    // Set the new particle position
	    particles[ppos_idx] = bpos;
	    // Remove the particle position from the border
	    std::swap(border[bpos_idx], border.back());
	    border.pop_back();
	    // Add the new borders
	    auto add_border = [&](long b) {
		if (grid[b] == -1 && b >= 2*(width+4) && b <= (height+2)*(width+4) && b%(width+4) >= 2 && b%(width+4) < width+2) {
		    bool add = true;
		    for (long pos : border) {
			if (pos == b) {
			    add = false;
			    break;
			}
		    }
		    if (add)
			border.push_back(b);
		}
	    };
	    add_border(ppos);
	    add_border(bpos-1);
	    add_border(bpos+1);
	    add_border(bpos-width-4);
	    add_border(bpos+width+4);
	    // return the energy difference
	    return dE;
	}
    }
}



void save_image(const long width, const long height, const char *grid, const std::vector<long>& border, const std::string& filename) {
    auto img = std::make_unique<char[]>((width+4)*(height+4));
    for (long i = 0 ; i < (width+4)*(height+4) ; ++i)
	img[i] = grid[i] == 1 ? 1 : 2;
    for (long b : border)
	img[b] = 0;
    std::ofstream imgstr(filename);
    imgstr << "P5\n" << width+4 << ' ' << height+4 << "\n2\n";
    for (long i = 0 ; i < (width+4)*(height+4) ; ++i)
	imgstr << img[i];
    imgstr.close();
}



int main(int argc, char* argv[]) {
    const double J = std::atof(argv[1]);
    const double K = std::atof(argv[2]);
    const double g = std::atof(argv[3]);
    const double T = std::atof(argv[4]);
    const long dropwidth = std::atol(argv[5]);
    const long dropheight = std::atol(argv[6]);
    const int reps = std::atol(argv[7]);
    const long width = 3*dropwidth/2;
    const long height = 3*dropheight;

    auto t = generate_grid(width, height, J, K, g, dropwidth, dropheight);
    auto grid = std::move(std::get<0>(t));
    auto particles = std::get<1>(t);
    auto border = std::get<2>(t);
    auto E = std::get<3>(t);

    std::mt19937 rng(0);
    std::cerr << std::string(reps, '-') << std::endl;
    for (int i = 0 ; i < reps ; ++i) {
	std::cerr << '-';
	char filename[30];
	sprintf(filename, "/tmp/anim/frame%04d.pgm", i);
	save_image(width, height, grid.get(), border, filename);
	for (int j = 0 ; j < 2000 ; ++j)
	    E += step(width, height, J, K, g, T, grid.get(), particles, border, rng);
    }
    std::cerr << std::endl;

    return 0;
}
