#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// Function to generate the polynomial values
double polynomial(double x) {
    return 0 + x + 2 * pow(x, 2) + 3 * pow(x, 3) + 4 * pow(x, 4) +
           5 * pow(x, 5) + 6 * pow(x, 6) + 7 * pow(x, 7);
}

// Function to generate matrix A and vector y
void generateMatrix(vector<vector<double>> &A, vector<vector<double>> &y, const vector<double> &x) {
    int rows = x.size();
    int cols = 8; // Degree 7 polynomial
    for (int i = 0; i < rows; ++i) {
        vector<double> row;
        double val = 1.0; // Start with x^0 = 1
        for (int j = 0; j < cols; ++j) {
            row.push_back(val);
            val *= x[i]; // Compute x^j
        }
        A.push_back(row);
        y.push_back({polynomial(x[i])}); // y is now a 2D matrix with one column
    }
}

// Function to multiply two matrices or a matrix and a vector
vector<vector<double>> multiply(const vector<vector<double>> &A, const vector<vector<double>> &B) {
    int rows = A.size(), cols = B[0].size(), common = A[0].size();
    vector<vector<double>> result(rows, vector<double>(cols, 0.0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < common; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Function to transpose a matrix
vector<vector<double>> transposeMatrix(const vector<vector<double>> &A) {
    int rows = A.size(), cols = A[0].size();
    vector<vector<double>> result(cols, vector<double>(rows));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = A[i][j];
        }
    }
    return result;
}

// Gaussian Elimination Solver
void gaussianElimination(vector<vector<double>> &B, vector<vector<double>> &d, vector<double> &c) {
    int n = B.size();
    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(B[k][i]) > abs(B[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(B[i], B[maxRow]);
        swap(d[i], d[maxRow]);

        // Elimination
        for (int k = i + 1; k < n; ++k) {
            double factor = B[k][i] / B[i][i];
            for (int j = i; j < n; ++j) {
                B[k][j] -= factor * B[i][j];
            }
            d[k][0] -= factor * d[i][0];
        }
    }

    // Back substitution
    c.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        c[i] = d[i][0];
        for (int j = i + 1; j < n; ++j) {
            c[i] -= B[i][j] * c[j];
        }
        c[i] /= B[i][i];
    }
}

int main() {
    // Step 1: Generate sample data
    vector<double> x = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    vector<vector<double>> A;
    vector<vector<double>> y;
    generateMatrix(A, y, x);

    // Step 2: Compute A^T * A and A^T * y
    vector<vector<double>> A_T = transposeMatrix(A);
    vector<vector<double>> B = multiply(A_T, A);
    vector<vector<double>> d = multiply(A_T, y);

    // Step 3: Solve Bc = d using Gaussian elimination
    vector<double> c;
    gaussianElimination(B, d, c);


    for(auto i:B) {
        for(auto j:i){
            cout << fixed << setprecision(10) << j << " ";
        }
        cout << endl;
    }
    // Output the solution
    cout << fixed << setprecision(10);
    cout << "Solution (coefficients of the polynomial):" << endl;
    for (size_t i = 0; i < c.size(); ++i) {
        cout << "a" << i << " = " << c[i] << endl;
    }

    return 0;
}
