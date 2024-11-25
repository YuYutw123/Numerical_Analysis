// QR Householder new

#include <bits/stdc++.h>
#define double long double
using namespace std;

struct poly {
    vector<double> coe;
};

double vector_norm(vector<double> &v) {
    double sum = 0;
    for(double val: v) {
        sum += val*val;
    }
    return sqrt(sum);
}

double infinity_norm(vector<double> &v) {
    double max_val = 0;
    for (double val : v) {
        max_val = max(max_val, fabs(val));
    }
    return max_val;
}

double function_p(poly p1, double x) {
    double sum = 0;
    for(int i=p1.coe.size()-1;i>=0;i--){
        sum = sum*x+p1.coe[i]; // Horner's method
    }
    return sum;
}

vector<vector<double>> build_matrix(vector<double> xi, int degree) {
    vector<vector<double>> mt(xi.size(), vector<double>(degree + 1, 1.0));
    for(int i=0;i<xi.size();i++) {
        for(int j=1;j<=degree;j++) {
            mt[i][j] = mt[i][j-1] *xi[i];
        }
    } 
    return mt;
}

vector<vector<double>> transpose_matrix(vector<vector<double>> &mt) {
    int n = mt.size(), m = mt[0].size();
    vector<vector<double>> tmp(m, vector<double>(n, 1.0));
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            tmp[j][i] = mt[i][j];
        }
    }
    return tmp;
}

vector<vector<double>> multiply_matrix(vector<vector<double>> A, vector<vector<double>> B) {
    int n = A.size(), m = B[0].size(), p = A[0].size();
    vector<vector<double>> tmp(n, vector<double>(m, 0.0));
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            for (int k=0;k<p;k++) {
                tmp[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return tmp;
}


void householder_qr(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size(), m = A[0].size();
    vector<vector<double>> H = A;

    Q = vector<vector<double>>(n, vector<double>(n, 0.0)); // Identity matrix
    for (int i = 0; i < n; ++i) {
        Q[i][i] = 1.0;
    }

    for (int k = 0; k < m; ++k) {
        // vector x
        vector<double> x(n - k, 0.0);
        for (int i = k; i < n; ++i) {
            x[i - k] = H[i][k];
        }

        // v and ||x||
        double norm_x = vector_norm(x);
        if (norm_x == 0.0) continue;
        x[0] += (x[0] >= 0 ? norm_x : -norm_x);
        double norm_v = vector_norm(x);
        for (double& val : x) {
            val /= norm_v;
        }

        // implement H_k
        vector<vector<double>> v_mat(n - k, vector<double>(1, 0.0));
        for (int i = 0; i < n - k; ++i) {
            v_mat[i][0] = x[i];
        }
        vector<vector<double>> v_t = transpose_matrix(v_mat);
        vector<vector<double>> v_vt = multiply_matrix(v_mat, v_t);
        vector<vector<double>> H_k = vector<vector<double>>(n - k, vector<double>(n - k, 0.0));
        for (int i = 0; i < n - k; ++i) {
            for (int j = 0; j < n - k; ++j) {
                H_k[i][j] = (i == j ? 1.0 : 0.0) - 2.0 * v_vt[i][j];
            }
        }

        // update H and Q
        vector<vector<double>> H_k_full(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i >= k && j >= k) {
                    H_k_full[i][j] = H_k[i - k][j - k];
                } else {
                    if (i == j) {
                        H_k_full[i][j] = 1.0;
                    } else {
                        H_k_full[i][j] = 0.0;
                    }
                }

            }
        }

        Q = multiply_matrix(Q, H_k_full);
        H = multiply_matrix(H_k_full, H);
    }
    R = H;
}
int main() {
    vector<double> xi = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    vector<vector<double>> y;
    vector<vector<double>> mt_Q, mt_R;
    poly p;
    for(int i=0;i<8;i++){
        p.coe.push_back(i);
    }
    for(int i=0;i<15;i++){
        // xi.push_back(val);
        y.push_back({function_p(p, xi[i])});
    }


    // AT * A * c = AT * y
    // B * c = d 
    // B = Q * R = AT * A
    // d = AT * y
    // Q * R * c = d
    // R * c = QT * d
    

    vector<vector<double>> mt = build_matrix(xi, 7);

    // AT * A => B
    vector<vector<double>> mt_T = transpose_matrix(mt);
    vector<vector<double>> mt_B = multiply_matrix(mt_T, mt);


    // AT * y => d
    vector<vector<double>> mt_d = multiply_matrix(mt_T, y);
    householder_qr(mt_B, mt_Q, mt_R);

    vector<vector<double>> Q_T = transpose_matrix(mt_Q);
    vector<vector<double>> right = multiply_matrix(Q_T, mt_d);

    vector<double> c(8);
    for (int i = mt_R.size() - 1; i >= 0; --i) {
        c[i] = right[i][0];
        for (int j = i + 1; j < mt_R[0].size(); ++j) {
            c[i] -= mt_R[i][j] * c[j];
        }
        c[i] /= mt_R[i][i];
    }
    for(auto i:c) {
        cout << fixed << setprecision(10) << i << endl;
    }
    cout << endl;
}