#include <bits/stdc++.h>
using namespace std;

struct poly {
    vector<double> coe;
};

// Horner's method for evaluating polynomial
double function_p(poly p1, double x) {
    double sum = 0;
    for (int i = p1.coe.size() - 1; i >= 0; i--) {
        sum = sum * x + p1.coe[i];
    }
    return sum;
}

// Build matrix for A
vector<vector<double>> build_matrix(const vector<double> &xi, int degree) {
    vector<vector<double>> mt(xi.size(), vector<double>(degree + 1, 1.0));
    for (int i = 0; i < xi.size(); i++) {
        for (int j = 1; j <= degree; j++) {
            mt[i][j] = mt[i][j - 1] * xi[i];
        }
    }
    return mt;
}

// Transpose matrix
vector<vector<double>> transpose_matrix(const vector<vector<double>> &mt) {
    int n = mt.size(), m = mt[0].size();
    vector<vector<double>> tmp(m, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            tmp[j][i] = mt[i][j];
        }
    }
    return tmp;
}

// Multiply two matrices
vector<vector<double>> multiply_matrix(const vector<vector<double>> &A, const vector<vector<double>> &B) {
    int n = A.size(), m = B[0].size(), p = A[0].size();
    vector<vector<double>> tmp(n, vector<double>(m, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                tmp[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return tmp;
}

// Gaussian elimination solver
vector<double> gauss_eliminate(vector<vector<double>> A, vector<vector<double>> B) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        // Partial pivoting
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[i][i])) {
                swap(A[j], A[i]);
                swap(B[j][0], B[i][0]);
            }
        }
        // Elimination
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            B[j][0] -= factor * B[i][0];
        }
    }
    // Back substitution
    vector<double> ans(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = B[i][0];
        for (int j = i + 1; j < n; j++) {
            ans[i] -= A[i][j] * ans[j];
        }
        ans[i] /= A[i][i];
    }
    return ans;
}

int main() {
    vector<double> xi, yi;
    poly p;
    for (int i = 0; i < 8; i++) {
        p.coe.push_back(i); // Polynomial coefficients
    }

    // Generate xi and yi
    double val = 0.2;
    for (int i = 0; i < 15; i++) {
        xi.push_back(val);
        yi.push_back(function_p(p, val));
        val += 0.2;
    }

    // Build matrix A
    vector<vector<double>> A = build_matrix(xi, 7);

    // Compute B = A^T * A
    vector<vector<double>> A_T = transpose_matrix(A);
    vector<vector<double>> B = multiply_matrix(A_T, A);

    // Compute d = A^T * y
    vector<vector<double>> y_vec(yi.size(), vector<double>(1, 0.0));
    for (int i = 0; i < yi.size(); i++) {
        y_vec[i][0] = yi[i];
    }
    vector<vector<double>> d = multiply_matrix(A_T, y_vec);

    // Solve for coefficients
    vector<double> c = gauss_eliminate(B, d);

    // Output the coefficients
    cout << "Polynomial coefficients:" << endl;
    for (auto coef : c) {
        cout << fixed << setprecision(10) << coef << " ";
    }
    cout << endl;

    return 0;
}
