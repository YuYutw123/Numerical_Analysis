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
    for(int i=0;i<15;i++) {
        cout << "(x" << i << ", y" << i << "): ";
        cout << "(" << xi[i] << ", " << yi[i] << ")" << endl;
    }
}