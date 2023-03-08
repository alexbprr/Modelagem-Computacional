#include <iostream>
#include <vector>
#include <random>

// Function to calculate propensity values for each reaction
std::vector<double> calc_propensity(double x, double y, double k1, double k2) {
    std::vector<double> a(2);
    a[0] = k1 * x;
    a[1] = k2 * x * y;
    return a;
}

// Function to update the system state based on a chosen reaction
void update_state(int& x, int& y, int reaction) {
    switch(reaction) {
        case 0: // Reaction 1
            x++;
            break;
        case 1: // Reaction 2
            x--;
            y--;
            break;
    }
}

// Function to choose the next reaction to occur based on the current propensity values
int choose_reaction(double a0, const std::vector<double>& a, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double r1 = dis(gen);
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += a[i];
        if (r1 < sum/a0) {
            return i;
        }
    }
    return -1; // No reaction chosen
}

// Main function to run the Gillespie algorithm
int main() {
    // Initialize system state and parameters
    int x = 2;
    int y = 8;
    double k1 = 0.1;
    double k2 = 0.05;
    double t = 0.0;
    int num_steps = 1000;
    
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> plot_time(num_steps);
    std::vector<int> plot_x(num_steps);
    std::vector<int> plot_y(num_steps);
    
    // Main loop
    for (int i = 0; i < num_steps; i++) {
        // Calculate propensity values for each reaction
        std::vector<double> a = calc_propensity(x, y, k1, k2);
        
        // Calculate the sum of the propensity values
        double a0 = std::accumulate(a.begin(), a.end(), 0.0);
        
        // Choose the next reaction to occur based on the current propensity values
        int reaction = choose_reaction(a0, a, gen);
        
        // If no reaction chosen, break out of loop
        if (reaction == -1) {
            std::cout << "No reaction chosen!" << std::endl;
            break;
        }
        
        // Update the system state based on the chosen reaction
        update_state(x, y, reaction);
        
        // Calculate the time elapsed for this step
        double r2 = dis(gen);
        double tau = 1.0/a0 * std::log(1.0/r2);
        
        // Update the time
        t += tau;
        
        // Output the current state and time
        std::cout << "Step " << i+1 << ": x = " << x << ", y = " << y << ", t = " << t << std::endl;

        plot_time.push_back(t);
        plot_x.push_back(x);
        plot_y.push_back(y);
    }

    // Save plot data to CSV file
    std::ofstream file("data.csv");
    if (file.is_open()) {
        file << "Time,X,Y" << std::endl;
        for (size_t i = 0; i < plot_x.size(); i++) {
            file << plot_time[i] << "," << plot_x[i] << "," << plot_y[i] << std::endl;
        }
        file.close();
        std::cout << "Data saved to data.csv" << std::endl;
    } else {
        std::cout << "Unable to open file" << std::endl;
    }
    
    return 0;
}