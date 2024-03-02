#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>

using namespace std;

// Finite difference coefficients
double alpha(int n, double sigma, double r, double D, double dt) {
    return 0.5 * (n * n * sigma * sigma - n * (r - D)) * dt;
}

double beta(int n, double sigma, double r, double D, double dt) {
    return 1 - (r + n * n * sigma * sigma) * dt;
}

double gamma(int n, double sigma, double r, double D, double dt) {
    return 0.5 * (n * n * sigma * sigma + n * (r - D)) * dt;
}

int main() {
    int S0 = 100; // Initial stock price
    int Smax = 3*S0; // Maximum stock price (Boundary condition)
    int E = 100; // Strike
    float Sigma = 0.2; // Volatility
    float r = 0.05; // Risk-free rate
    float D = 0.0; // Dividend yield
    int T = 1; // Time to maturity
    int tgrid = 100000; // No. time grid points
    int sgrid = 1000; // No. price grid points
    double dt = 1.0 * T/tgrid;
    double ds = 1.0 * Smax/sgrid;
    int index = S0/ds; // Closest index to initial price
    const char* type = "American";
    const char* direction = "Put";

    // Check for stability
    if (dt >= 1/(Sigma*Sigma * sgrid*sgrid)){
        cout << "Unstable!" << endl;
        exit(0);
    }
    
    // Interpolate across grid for intial stock price
    float leftdiff = abs(ds * index - S0);
    float rightdiff = abs(ds * (index+1) - S0);
    float leftweight = rightdiff/(leftdiff+rightdiff);
    float rightweight = leftdiff/(leftdiff+rightdiff);
    
    // Intialise grid
    const int rows = sgrid;
    const int cols = tgrid;
    vector<vector<float>> grid(rows, vector<float>(cols, 0.0)); // Index by: grid[row][column] or grid[sgrid][tgrid]
    
    

    
    
    
    // Call Option Initial Values
    if (strcmp(direction, "Call") == 0){
        
        // Final Payoff Boundary conditions
        for (int i = 0; i < sgrid; i++){
            
            float val = std::max<float>(i*ds-E, 0);
            grid[i][tgrid-1] = val;
        }
        
        // S = 0 Boundary Conditions
        for (int i = 0; i < tgrid; i++){
            
            grid[0][i] = 0;
            
        }
        
        // S = Smax Boundary Conditions
        for (int i = 0; i < tgrid; i++){
            
            grid[sgrid-1][i] = max(Smax - E - ds, 0.0);
            
        }
        
    }
    
    // Put Option Initial Values
    if (strcmp(direction, "Put") == 0){
        
        // Final Payoff Boundary conditions
        for (int i = 0; i < sgrid; i++){
            
            float val = std::max<float>(E - i*ds, 0);
            grid[i][tgrid-1] = val;
        }
        
        // S = 0 Boundary Conditions
        for (int i = 0; i < tgrid; i++){
            
            grid[0][i] = E;
            
        }
        
        // S = Smax Boundary Conditions
        for (int i = 0; i < tgrid; i++){
            
            grid[sgrid-1][i] = 0;
            
        }
        
    }
    

    // Finite Difference Updates:
    
    // Backwards Marching Time Loop
    for (int m = tgrid - 2; m >= 0; m--) {
        
        // Stock price loop
        for (int n = 1; n < sgrid - 1; n++) {
            
            grid[n][m] = alpha(n, Sigma, r, D, dt) * grid[n - 1][m + 1] +
            beta(n, Sigma, r, D, dt) * grid[n][m + 1] +
            gamma(n, Sigma, r, D, dt) * grid[n + 1][m + 1];
        }
        
        // Boundary Condition at N (With Gamma/V" = 0)
        grid[sgrid-1][m] = (alpha(sgrid-1, Sigma, r, D, dt) - gamma(sgrid-1, Sigma, r, D, dt)) * grid[sgrid-2][m+1] + (beta(sgrid-1, Sigma, r, D, dt) + 2 * gamma(sgrid-1, Sigma, r, D, dt)) * grid[sgrid-1][m+1];
        
        // American option condition
        if (strcmp(type, "American") == 0){
            for (int n = 0; n < sgrid-1; n++){
                if (grid[n][m] < grid[n][tgrid-1]){
                    grid[n][m] = grid[n][tgrid-1];
                }
            }
        }
        
    }
    
    // Interpolation
    float result = leftweight * grid[index][0] + rightweight * grid[index+1][0];
    
    cout << "Option Price: " << result << "\n" << endl;
    
    return 0;
}
