
#include "HEaaN/HEaaN.hpp"
#include "ccmm.hpp"
#include "math.h"
#include <cstdint>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <iostream>
#include "HEaaNTimer.hpp"
#include "HEaaN-math/HEaaN-math.hpp"
#include "examples.hpp"
#include <map>
#include <algorithm>


int main() {

    int d = 128;
    int dd = d*d;
    int slot = 2*dd;

    std::cout << "d : " << d << ", dd : " << dd << ", slot : " << slot << '\n';
   
    std::unordered_map<int, std::vector<double>> sigma,tau,transpose;
    std::vector<std::unordered_map<int, std::vector<double>>> phi;

    if (!std::filesystem::exists("sigma.bin")) {
        sigma = make_sigma(d,dd,slot);

        try {
            saveBinary("sigma.bin", sigma);
            std::cout << "바이너리 저장 완료\n";

            auto loaded = loadBinary("sigma.bin");
            std::cout << "로드한 데이터 크기: " << loaded.size() << "\n";
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
        }
    }

    sigma = loadBinary("sigma.bin");

    for(int i=0; i<10; i++) {
        std::cout << sigma[0][i] << " ";
    }
    std::cout << '\n';

    if (!std::filesystem::exists("tau.bin")) {
        tau = make_tau(d,dd,slot);

        try {
            saveBinary("tau.bin", tau);
            std::cout << "바이너리 저장 완료\n";

            auto loaded = loadBinary("tau.bin");
            std::cout << "로드한 데이터 크기: " << loaded.size() << "\n";
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
        }
    }

    tau = loadBinary("tau.bin");

    for(int i=0; i<10; i++) {
        std::cout << tau[0][i] << " ";
    }
    std::cout << '\n';

    if (!std::filesystem::exists("transpose.bin")) {
        transpose = make_transpose(d,dd,slot);

        try {
            saveBinary("transpose.bin", transpose);
            std::cout << "바이너리 저장 완료\n";

            auto loaded = loadBinary("transpose.bin");
            std::cout << "로드한 데이터 크기: " << loaded.size() << "\n";
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
        }
    }

    transpose = loadBinary("transpose.bin");

    for(int i=0; i<10; i++) {
        std::cout << transpose[0][i] << " ";
    }
    std::cout << '\n';

    if (!std::filesystem::exists("phi.bin")) {
        phi = make_phi(d,dd,slot);

        try {
            saveBinary2("phi.bin", phi);
            std::cout << "바이너리 저장 완료\n";

            auto loaded = loadBinary2("phi.bin");
            std::cout << "로드한 데이터 크기: " << loaded.size() << "\n";
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
        }
    }

    phi = loadBinary2("phi.bin");

    std::cout << "go\n";
    for(int i=0; i<10; i++) {
        std::cout << phi[0][0][i] << " ";
    }
    std::cout << '\n';


    std::unordered_map<int, std::vector<std::complex<double>>> ext_sigma;
    for (const std::pair<const int, std::vector<double>>& p : sigma) {
        int key   = p.first;
        std::vector<double> val = p.second;

        std::vector<std::complex<double>> temp;
        for(int i=0; i<val.size(); i++) {
            temp.push_back(std::complex<double>(val[i],0));
        }
        ext_sigma[key] = temp;
    }

    saveBinaryComplex("complex_sigma.bin",ext_sigma);

    std::unordered_map<int, std::vector<std::complex<double>>> ext_tau;
    for (const std::pair<const int, std::vector<double>>& p : tau) {
        int key   = p.first;
        std::vector<double> val = p.second;

        std::vector<std::complex<double>> temp;
        for(int i=0; i<val.size(); i++) {
            temp.push_back(std::complex<double>(val[i],0));
        }
        ext_tau[key] = temp;
    }

    saveBinaryComplex("complex_tau.bin",ext_tau);
    
    int count = 0;
    std::vector<int> v;
    std::unordered_map<int, std::vector<std::complex<double>>> ext_transpose;
    for (const std::pair<const int, std::vector<double>>& p : transpose) {
        int key   = p.first;
        std::vector<double> val = p.second;
        count ++;
        //std::cout << key << '\n';
        v.push_back(key);
        std::vector<std::complex<double>> temp;
        for(int i=0; i<val.size(); i++) {
            temp.push_back(std::complex<double>(val[i],0));
        }
        ext_transpose[key] = temp;
    }
    std::cout << "count : " << count << '\n';
    sort(v.begin(),v.end());
    std::cout << "size: "  << v.size() << '\n';
    for(int i=0; i<v.size(); i++) {
        std::cout << v[i] << '\n';
    }
    saveBinaryComplex("complex_transpose.bin",ext_transpose);

    std::vector<std::unordered_map<int, std::vector<std::complex<double>>>> ext_phi(phi.size());
    
    for(int k=0; k<phi.size(); k++) {
        
        for (const std::pair<const int, std::vector<double>>& p : phi[k]) {
            int key  = p.first;
            std::vector<double> val = p.second;
    
            std::vector<std::complex<double>> temp;
            for(int i=0; i<val.size(); i++) {
                temp.push_back(std::complex<double>(val[i],0));
            }
            ext_phi[k][key] = temp;
        }
    }

    saveBinary2Complex("complex_phi.bin",ext_phi);

}