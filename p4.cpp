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

int main() {
    vector<double> xi, yi;
    poly p;
    for(int i=0;i<8;i++){
        p.coe.push_back(i);
    }
    double val = 0.2;
    for(int i=0;i<15;i++){
        xi.push_back(val);
        yi.push_back(function_p(p, val));
        val += 0.2;
    }
    vector<vector<double>> mt = build_matrix(xi, 7);
    for(auto i: mt) {
        for (auto j: i) {
            cout << j << " ";
        }
        cout << endl;
    }
}