// Gauss Eliminate with partial pivoting

#include <bits/stdc++.h>
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
    vector<double> xi = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    vector<double> ori_c = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    int degree = 7;
    vector<vector<double>> result;
    vector<vector<double>> errors;
    vector<double> norms_2, norms_inf;
    ofstream out("data.txt");
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

        vector<double> error(mt_c.size());
        for (int i = 0; i < mt_c.size(); ++i) {
            error[i] = mt_c[i] - ori_c[i];
        }

        double norm_2 = vector_norm(error);
        double norm_inf = infinity_norm(error);

        errors.push_back(error);
        norms_2.push_back(norm_2);
        norms_inf.push_back(norm_inf);
    }
    out << "\nComparison of Norms:\n";
    cout << "\nComparison of Norms:\n";
    out << "Degree\t2-Norm\t\tInfinity-Norm\n";
    cout << "Degree\t2-Norm\t\tInfinity-Norm\n";
    for (int i = 0; i < norms_2.size(); ++i) {
        out << fixed << setprecision(15) << (5 + i) << "\t" << norms_2[i] << "\t" << norms_inf[i] << endl;
        cout << fixed << setprecision(15) << (5 + i) << "\t" << norms_2[i] << "\t" << norms_inf[i] << endl;
    }
    out.close();
}