#include<bits/stdc++.h>
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

vector<vector<double>> transpose_matrix(vector<vector<double>> &mt, int n, int m) {
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

int main() {
    vector<double> xi;
    vector<vector<double>> y;
    poly p;
    for(int i=0;i<8;i++){
        p.coe.push_back(i);
    }
    double val = 0.2;
    for(int i=0;i<15;i++){
        xi.push_back(val);
        y.push_back({function_p(p, val)});
        val += 0.2;
    }
    
    // AT * A * c = AT * y
    // B * c = d

    vector<vector<double>> mt = build_matrix(xi, 7);

    // AT * A => B
    vector<vector<double>> mt_T = transpose_matrix(mt,15, 8);
    vector<vector<double>> mt_B = multiply_matrix(mt_T, mt);


    // AT * y => d
    vector<vector<double>> mt_d = multiply_matrix(mt_T, y);
    vector<double> mt_c;
    gauss_eliminate(mt_B, mt_d, mt_c);
    
    for(auto i:mt_B) {
        for(auto j:i){
            cout << fixed << setprecision(10) << j << " ";
        }
        cout << endl;
    }
    for(auto i: mt_c) {
        cout << i << " ";
    }
    cout << endl;
}