// Gauss Eliminate with partial pivoting

#include <bits/stdc++.h>
#define double float
using namespace std;

struct poly {
    vector<double> coe;
};

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

void gauss_eliminate(vector<vector<double>> &A,vector<vector<double>> &B, vector<double> &c) {
    for(int i=0;i<A.size();i++) {
        for(int j=i+1;j<A.size();j++) {
            if(fabs(A[j][i]) > fabs(A[i][i])) {
                swap(A[j], A[i]);
                swap(B[j], B[i]);
            }
        }
        for(int j=i+1;j<A.size();j++) {
            double tmp = A[j][i]/A[i][i];
            for(int k=i;k<A.size();k++) {
                A[j][k] -= tmp*A[i][k];
            }
            B[j][0] -= tmp*B[i][0];
        }
    }
    c.resize(A.size());
    for(int i=A.size()-1;i>=0;i--) {
        c[i] = B[i][0];
        for(int j=i+1;j<A.size();j++) {
            c[i] -= A[i][j]*c[j];
        }
        c[i] /= A[i][i];
    }
}

double relative_error(const vector<double>& computed, const vector<double>& original) {
    double error = 0.0;
    for (size_t i = 0; i < computed.size(); ++i) {
        error += (computed[i] - original[i]) * (computed[i] - original[i]);
    }
    error = sqrt(error);

    double norm_original = 0.0;
    for (double val : original) {
        norm_original += val * val;
    }
    norm_original = sqrt(norm_original);

    return error / norm_original;
}

int main() {
    vector<double> xi = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    vector<double> ori_c = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    int degree = 7;
    vector<vector<double>> result;
    vector<double> rel_error;
    for(degree=5;degree<=8;degree++) {
        vector<vector<double>> y;
        vector<vector<double>> mt_Q, mt_R;
        poly p;
        for(int i=0;i<degree+1;i++){
            p.coe.push_back(i);
        }

        for(int i=0;i<15;i++){
            // xi.push_back(val);
            y.push_back({function_p(p, xi[i])});
        }
        
        // AT * A * c = AT * y
        // B * c = d

        vector<vector<double>> mt = build_matrix(xi, degree);

        // AT * A => B
        vector<vector<double>> mt_T = transpose_matrix(mt);
        vector<vector<double>> mt_B = multiply_matrix(mt_T, mt);


        // AT * y => d
        vector<vector<double>> mt_d = multiply_matrix(mt_T, y);
        vector<double> mt_c;
        gauss_eliminate(mt_B, mt_d, mt_c);
        result.push_back(mt_c);
        rel_error.push_back(relative_error(mt_c, ori_c));
        
        // Check coeff
        // for(auto i: mt_c) {
        //     cout << fixed << setprecision(20) << i << " ";
        // }
        // cout << endl;
    }

    int cnt = 5;
    for(auto i: result) {
        cout << "Degree " << cnt << ": ";
        for(auto j:i){
            cout << fixed << setprecision(10) << j << " ";
        }
        cout << endl;
        cout << "Relative Error: ";
        cout << rel_error[cnt-5] << endl;
        cnt++;
    }
    
}