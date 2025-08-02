#pragma once

#include "HEaaN/HEaaN.hpp"
#include "HEaaN-math/HEaaN-math.hpp"
#include <unordered_map>
#include "examples.hpp"
#include <algorithm>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cstdint>
#include <string>

void saveBinaryComplex(
    const std::string& filename,
    const std::unordered_map<int, std::vector<std::complex<double>>>& m)
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("파일 열기 실패: " + filename);
    }

    // 1) 맵 크기 저장
    uint64_t mapSize = m.size();
    ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

    // 2) 각 엔트리 순회
    for (const auto& entry : m) {
        int key = entry.first;
        const std::vector<std::complex<double>>& vec = entry.second;

        // key 저장
        ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
        // 벡터 길이 저장
        uint64_t vecSize = vec.size();
        ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));
        // complex<double> 데이터 저장
        ofs.write(
            reinterpret_cast<const char*>(vec.data()),
            vecSize * sizeof(std::complex<double>)
        );
    }

    ofs.close();
}

std::unordered_map<int, std::vector<std::complex<double>>>
loadBinaryComplex(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("파일 열기 실패: " + filename);
    }

    // 1) 맵 크기 읽기
    uint64_t mapSize;
    ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

    std::unordered_map<int, std::vector<std::complex<double>>> m;
    for (uint64_t i = 0; i < mapSize; ++i) {
        // key 읽기
        int key;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));

        // 벡터 길이 읽기
        uint64_t vecSize;
        ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));

        // complex<double> 벡터 읽기
        std::vector<std::complex<double>> vec(vecSize);
        ifs.read(
            reinterpret_cast<char*>(vec.data()),
            vecSize * sizeof(std::complex<double>)
        );

        m.emplace(key, std::move(vec));
    }

    ifs.close();
    return m;
}


void saveBinary(const std::string& filename,
    const std::unordered_map<int, std::vector<double>>& m) {
std::ofstream ofs(filename, std::ios::binary);
if (!ofs) {
throw std::runtime_error("파일 열기 실패: " + filename);
}

// 1) 맵 크기 저장
uint64_t mapSize = m.size();
ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

// 2) 각 엔트리 순회
for (auto const& [key, vec] : m) {
// key 저장
ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
// 벡터 길이 저장
uint64_t vecSize = vec.size();
ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));
// 벡터 데이터 저장
ofs.write(reinterpret_cast<const char*>(vec.data()),
      vecSize * sizeof(double));
}
ofs.close();
}

std::unordered_map<int, std::vector<double>>
loadBinary(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("파일 열기 실패: " + filename);
    }

    // 1) 맵 크기 읽기
    uint64_t mapSize;
    ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

    std::unordered_map<int, std::vector<double>> m;
    for (uint64_t i = 0; i < mapSize; ++i) {
        int key;
        ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
        uint64_t vecSize;
        ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
        std::vector<double> vec(vecSize);
        ifs.read(reinterpret_cast<char*>(vec.data()),
                 vecSize * sizeof(double));
        m.emplace(key, std::move(vec));
    }
    ifs.close();
    return m;
}

// 저장
void saveBinary2(
    const std::string& filename,
    const std::vector<std::unordered_map<int, std::vector<double>>>& vmap)
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("파일 열기 실패: " + filename);

    // 벡터 크기
    uint64_t vecSize = vmap.size();
    ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));

    for (size_t i = 0; i < vmap.size(); ++i) {
        const std::unordered_map<int, std::vector<double>>& mp = vmap[i];
        uint64_t mapSize = mp.size();
        ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

        for (std::unordered_map<int, std::vector<double>>::const_iterator it = mp.begin();
             it != mp.end(); ++it)
        {
            int key = it->first;
            const std::vector<double>& vec = it->second;
            // key
            ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
            // 벡터 길이
            uint64_t vlen = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vlen), sizeof(vlen));
            // 벡터 데이터
            ofs.write(reinterpret_cast<const char*>(vec.data()),
                      vlen * sizeof(double));
        }
    }
}

// 로드
std::vector<std::unordered_map<int, std::vector<double>>>
loadBinary2(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("파일 열기 실패: " + filename);

    uint64_t vecSize;
    ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));

    std::vector<std::unordered_map<int, std::vector<double>>> vmap;
    vmap.reserve(vecSize);

    for (uint64_t i = 0; i < vecSize; ++i) {
        uint64_t mapSize;
        ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

        std::unordered_map<int, std::vector<double>> mp;
        for (uint64_t k = 0; k < mapSize; ++k) {
            int key;
            ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
            uint64_t vlen;
            ifs.read(reinterpret_cast<char*>(&vlen), sizeof(vlen));
            std::vector<double> vec(vlen);
            for (size_t j = 0; j < vec.size(); ++j) {
                ifs.read(reinterpret_cast<char*>(&vec[j]), sizeof(double));
            }
            mp.emplace(key, std::move(vec));
        }
        vmap.push_back(std::move(mp));
    }
    return vmap;
}

void saveBinary2Complex(
    const std::string& filename,
    const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& vmap)
{
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("파일 열기 실패: " + filename);

    // 벡터 크기
    uint64_t vecSize = vmap.size();
    ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));

    for (size_t i = 0; i < vmap.size(); ++i) {
        const std::unordered_map<int, std::vector<std::complex<double>>>& mp = vmap[i];
        uint64_t mapSize = mp.size();
        ofs.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));

        for (std::unordered_map<int, std::vector<std::complex<double>>>::const_iterator it = mp.begin();
             it != mp.end(); ++it)
        {
            int key = it->first;
            const std::vector<std::complex<double>>& vec = it->second;
            // key 저장
            ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));
            // 벡터 길이 저장
            uint64_t vlen = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vlen), sizeof(vlen));
            // complex<double> 데이터 저장
            ofs.write(
                reinterpret_cast<const char*>(vec.data()),
                vlen * sizeof(std::complex<double>)
            );
        }
    }
}

std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>
loadBinary2Complex(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) throw std::runtime_error("파일 열기 실패: " + filename);

    uint64_t vecSize;
    ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));

    std::vector<std::unordered_map<int, std::vector<std::complex<double>>>> vmap;
    vmap.reserve(vecSize);

    for (uint64_t i = 0; i < vecSize; ++i) {
        uint64_t mapSize;
        ifs.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));

        std::unordered_map<int, std::vector<std::complex<double>>> mp;
        for (uint64_t k = 0; k < mapSize; ++k) {
            int key;
            ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
            uint64_t vlen;
            ifs.read(reinterpret_cast<char*>(&vlen), sizeof(vlen));
            std::vector<std::complex<double>> vec(vlen);
            ifs.read(
                reinterpret_cast<char*>(vec.data()),
                vlen * sizeof(std::complex<double>)
            );
            mp.emplace(key, std::move(vec));
        }
        vmap.push_back(std::move(mp));
    }
    return vmap;
}

std::vector<std::vector<double>> double_convert(const std::vector<std::vector<double>>& mat) {
    long len = mat.size();
    std::vector<std::vector<double>> rtn(2*len, std::vector<double>(2*len,0));

    for(int i=0; i<len; i++) {
        for(int j=0; j<len; j++) {
            rtn[2*i][2*j] = mat[i][j];
            rtn[2*i+1][2*j+1] = mat[i][j];
        }
    }

    return rtn;
}

std::unordered_map<int, std::vector<double>> make_sigma(int d, int dd, int slot) {
    // 1) d x d 매트릭스 생성

    std::cout << "debug\n";
    std::vector<std::vector<double>> mat(dd, std::vector<double>(dd,0.0));

    // 2) mat 채우기
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            int row = d * i + j;
            for (int L = 0; L < dd; ++L) {
                if (L == d * i + (i + j) % d) {
                    mat[row][L] = 1.0;
                } else {
                    mat[row][L] = 0.0;
                }
            }
        }
    }

    // 3) DoubleConvert 호출
    std::cout << "debug\n";
    std::vector<std::vector<double>> tmp = double_convert(mat);

    std::cout << "debug\n";
    // 4) 대각 성분만 골라 diagmat에 저장
    std::unordered_map<int, std::vector<double>> diagmat;
    for (int i = 0; i < tmp.size(); ++i) {
        std::vector<double> diag(slot,0);
        double sum = 0.0;
        for (int j=0; j<tmp.size(); ++j) {
            diag[j] = tmp[j][(j+i)%slot];
            sum += diag[j];
        }

        if(sum != 0.0) {
            diagmat[i] = diag;
        }
    }

    return diagmat;
}

std::unordered_map<int, std::vector<double>> make_tau(int d, int dd, int slot) {
    // 1) d x d 매트릭스 생성
    std::vector<std::vector<double>> mat(dd, std::vector<double>(dd));

    // 2) mat 채우기
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            int row = d * i + j;
            for (int L = 0; L < dd; ++L) {
                if (L == d*((i+j)%d)+j) {
                    mat[row][L] = 1.0;
                } else {
                    mat[row][L] = 0.0;
                }
            }
        }
    }

    // 3) DoubleConvert 호출
    std::vector<std::vector<double>> tmp = double_convert(mat);

    // 4) 대각 성분만 골라 diagmat에 저장
    std::unordered_map<int, std::vector<double>> diagmat;
    for (int i = 0; i < tmp.size(); ++i) {
        std::vector<double> diag(slot,0);
        double sum = 0.0;
        for (int j=0; j<tmp.size(); ++j) {
            diag[j] = tmp[j][(j+i)%slot];
            sum += diag[j];
        }

        if(sum != 0.0) {
            diagmat[i] = diag;
        }
    }

    return diagmat;
}

std::unordered_map<int, std::vector<double>> make_transpose(int d, int dd, int slot) {
    // 1) d x d 매트릭스 생성
    std::vector<std::vector<double>> mat(dd, std::vector<double>(dd));

    // 2) mat 채우기
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            int row = d * i + j;
            for (int L = 0; L < dd; ++L) {
                if (L == d*j+i) {
                    mat[row][L] = 1.0;
                } else {
                    mat[row][L] = 0.0;
                }
            }
        }
    }

    // 3) DoubleConvert 호출
    std::vector<std::vector<double>> tmp = double_convert(mat);

    // 4) 대각 성분만 골라 diagmat에 저장
    std::unordered_map<int, std::vector<double>> diagmat;
    for (int i = 0; i < tmp.size(); ++i) {
        std::vector<double> diag(slot,0);
        double sum = 0.0;
        for (int j=0; j<tmp.size(); ++j) {
            diag[j] = tmp[j][(j+i)%slot];
            sum += diag[j];
        }

        if(sum != 0.0) {
            diagmat[i] = diag;
        }
    }

    return diagmat;
}

std::vector<std::unordered_map<int, std::vector<double>>> make_phi(int d, int dd, int slot) {
    
    std::vector<std::unordered_map<int, std::vector<double>>> rtn;

    for (int k=0; k<d; k++) {
        std::cout << "k : " << k << '\n';
         // 1) d x d 매트릭스 생성
        std::vector<std::vector<double>> mat(dd, std::vector<double>(dd));

        // 2) mat 채우기
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                int row = d * i + j;
                for (int L = 0; L < dd; ++L) {
                    if (L ==d*i+(j+k)%d ) {
                        mat[row][L] = 1.0;
                    } else {
                        mat[row][L] = 0.0;
                    }
                }
            }
        }

        // 3) DoubleConvert 호출
        std::vector<std::vector<double>> tmp = double_convert(mat);

        // 4) 대각 성분만 골라 diagmat에 저장
        std::unordered_map<int, std::vector<double>> diagmat;
        for (int i = 0; i < tmp.size(); ++i) {
            std::vector<double> diag(slot,0);
            double sum = 0.0;
            for (int j=0; j<tmp.size(); ++j) {
                diag[j] = tmp[j][(j+i)%slot];
                sum += diag[j];
            }

            if(sum != 0.0) {
                diagmat[i] = diag;
            }
        }

        rtn.push_back(diagmat);
    }
    
    return rtn;
}

std::vector<int> decomposeBinaryRemainder(unsigned int n) {
    std::vector<int> bits;
    if (n == 0) {
        bits.push_back(0);
        return bits;
    }
    while (n > 0) {
        bits.push_back(n % 2);  // LSB부터 차례대로
        n /= 2;
    }
    return bits;  // bits[0]은 2^0 자리, bits[1]은 2^1 자리, …
}

void leftRotate2(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt, int rot, HEaaN::Ciphertext &ctxt_res) {
    if (rot <0) {
        rot = (rot + 32768)%32768;
    }

    if (rot == 0 ){
        ctxt_res = ctxt;
        return;
    }
    
    //std::cout << rot << '\n';
    std::vector<int> decomposed = decomposeBinaryRemainder(rot); 

    
    bool start = true;
    HEaaN::Ciphertext temp(ctxt);
    for (int i=0; i<decomposed.size(); i++) {
        //std::cout << decomposed[i] << '\n';
        if (decomposed[i]==0) {
            continue;
        }
        //std::cout << (1<<i) << '\n';
        if(start) {
            eval.leftRotate(temp,1<<i,ctxt_res);
            start = false;
        } else {
            eval.leftRotate(ctxt_res,1<<i,ctxt_res);
        }
    }
}

#include <map>
#include "HEaaNTimer.hpp"


void sigma(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt, HEaaN::Ciphertext &ctxt_res, std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &plain) {

    //std::unordered_map<int, std::vector<std::complex<double>>> ext_sigma = loadBinaryComplex("complex_sigma.bin");
    int d = 128; // Hard coding
    int root2d = int(ceil(sqrt(2*d)));
    int slot = d*d;

    /*
    std::unordered_map<int,std::unique_ptr<HEaaN::Message>> plain;
    for(int j=0; j<root2d; j++) {
        for(int i=0; i<root2d; i++) {
            int shift = (-2*j + 2*slot)%(2*slot);
            std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);
           
            plain[2*((i*root2d+j-d+slot)%slot)] = std::make_unique<HEaaN::Message>();
            *plain[2*((i*root2d+j-d+slot)%slot)] = msg;
            plain[2*((i*root2d+j-d+slot)%slot)]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    std::cout << "Debug\n";
    */

    HEaaN::Ciphertext temp = HEaaN::Ciphertext(ctxt);

    // baby
    std::vector<uint64_t> giant_idx;
    for(int i=0; i< root2d; i++) {
        int index = (i*root2d-d+slot) % slot;
        giant_idx.push_back(2*index);
        //std::cout << 2*index << '\n';
    }

    //sort(giant_idx.begin(),giant_idx.end());

    //std::cout << "debug\n";

    std::unordered_map<int, std::unique_ptr<HEaaN::Ciphertext>> giant;

    // giant_idx가 std::vector<int>라 가정
    for (int i=0; i<giant_idx.size(); i++ ) {
        //std::cout << giant_idx[i]<< '\n';
        // 1) 새로운 Ciphertext 객체를 생성해 unique_ptr에 할당
        giant[giant_idx[i]] = std::make_unique<HEaaN::Ciphertext>(eval.getContext());
        giant[giant_idx[i]]->to(HEaaN::getCurrentCudaDevice());
        // 2) leftRotate2 호출 시, 포인터가 가리키는 객체를 참조로 전달
        if(i==0) 
            eval.leftRotate(temp, giant_idx[i], *giant[giant_idx[i]]);
        else
            eval.leftRotate(*giant[giant_idx[i-1]], 2*root2d, *giant[giant_idx[i]]);
            //leftRotate2(eval, *giant[giant_idx[i-1]], giant_idx[i]-giant_idx[i-1], *giant[giant_idx[i]]);
    }


    

    for(int j=0; j<root2d; j++) {
        HEaaN::Ciphertext inner(ctxt);
        inner.to(HEaaN::getCurrentCudaDevice());
        bool start = true;

        for(int i=0; i<root2d; i++) {

            if (!(-d<i*root2d+j-d && i*root2d+j-d<d)) {
                continue;
            }

            //int shift = (-2*j + 2*slot)%(2*slot);
            //std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            //std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            //HEaaN::Message msg(roted);
            //msg.to(HEaaN::getCurrentCudaDevice());

            int index = (i*root2d - d + slot) % slot;
            HEaaN::Ciphertext baby(temp);   
           

            eval.multWithoutRescale(*giant[2*index],*plain[2*((i*root2d+j-d+slot)%slot)],baby);
         
            if (start) {
                inner = baby;
                start = false;
                
            } else {
                eval.add(inner,baby,inner);
               
                //if(i==8) {
                //    ctxt_res = *giant[2*index];
                //    std::cout << 2*index << '\n';
                //    break;
               // }
                    
            }
        }

        //break;
        
        eval.rescale(inner);

      
        leftRotate2(eval,inner,2*j,inner);


        if (j==0) {
            ctxt_res = inner;
            //return;
        } else {
            eval.add(ctxt_res,inner,ctxt_res);
        }
    }

    
}


void tau(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt, HEaaN::Ciphertext &ctxt_res, std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &plain) {

    //std::unordered_map<int, std::vector<std::complex<double>>> ext_sigma = loadBinaryComplex("complex_sigma.bin");
    int d = 128; // Hard coding
    int rootd = int(ceil(sqrt(d)));
    int slot = d*d;

    /*
    std::unordered_map<int,std::unique_ptr<HEaaN::Message>> plain;
    for(int j=0; j<root2d; j++) {
        for(int i=0; i<root2d; i++) {
            int shift = (-2*j + 2*slot)%(2*slot);
            std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);
           
            plain[2*((i*root2d+j-d+slot)%slot)] = std::make_unique<HEaaN::Message>();
            *plain[2*((i*root2d+j-d+slot)%slot)] = msg;
            plain[2*((i*root2d+j-d+slot)%slot)]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    std::cout << "Debug\n";
    */

    HEaaN::Ciphertext temp = HEaaN::Ciphertext(ctxt);

    // baby
    std::vector<uint64_t> giant_idx;
    for(int i=0; i< rootd; i++) {
        int index = i*d*rootd;
        giant_idx.push_back(2*index);
        //std::cout << 2*index << '\n';
    }

    //sort(giant_idx.begin(),giant_idx.end());

    //std::cout << "debug\n";

    std::unordered_map<int, std::unique_ptr<HEaaN::Ciphertext>> giant;

    // giant_idx가 std::vector<int>라 가정
    for (int i=0; i<giant_idx.size(); i++ ) {
        //std::cout << giant_idx[i]<< '\n';
        // 1) 새로운 Ciphertext 객체를 생성해 unique_ptr에 할당
        giant[giant_idx[i]] = std::make_unique<HEaaN::Ciphertext>(eval.getContext());
        giant[giant_idx[i]]->to(HEaaN::getCurrentCudaDevice());
        // 2) leftRotate2 호출 시, 포인터가 가리키는 객체를 참조로 전달
        if(i==0) 
            eval.leftRotate(temp, giant_idx[i], *giant[giant_idx[i]]);
        else
            eval.leftRotate(*giant[giant_idx[i-1]], 2*d*rootd, *giant[giant_idx[i]]);
            //leftRotate2(eval, *giant[giant_idx[i-1]], giant_idx[i]-giant_idx[i-1], *giant[giant_idx[i]]);
    }


    

    for(int j=0; j<rootd; j++) {
        HEaaN::Ciphertext inner(ctxt);
        inner.to(HEaaN::getCurrentCudaDevice());
        bool start = true;

        for(int i=0; i<rootd; i++) {

            if (!(((i*rootd+j)*d)%d == 0 && (i*rootd+j)*d < slot)) {
                continue;
            }

            //int shift = (-2*j + 2*slot)%(2*slot);
            //std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            //std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            //HEaaN::Message msg(roted);
            //msg.to(HEaaN::getCurrentCudaDevice());

            int index = i * d * rootd;
            HEaaN::Ciphertext baby(temp);   
           

            eval.multWithoutRescale(*giant[2*index],*plain[2*(i*rootd+j)*d],baby);
         
            if (start) {
                inner = baby;
                start = false;
                
            } else {
                eval.add(inner,baby,inner);
               
                //if(i==8) {
                //    ctxt_res = *giant[2*index];
                //    std::cout << 2*index << '\n';
                //    break;
               // }
                    
            }
        }

        //break;
        
        eval.rescale(inner);

      
        leftRotate2(eval,inner,2*j*d,inner);


        if (j==0) {
            ctxt_res = inner;
            //return;
        } else {
            eval.add(ctxt_res,inner,ctxt_res);
        }
    }
}

void transpose(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt, HEaaN::Ciphertext &ctxt_res, std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &plain) {

    //std::unordered_map<int, std::vector<std::complex<double>>> ext_sigma = loadBinaryComplex("complex_sigma.bin");
    int d = 128; // Hard coding
    int root2d = int(ceil(sqrt(2*d)));
    int slot = d*d;

    /*
    std::unordered_map<int,std::unique_ptr<HEaaN::Message>> plain;
    for(int j=0; j<root2d; j++) {
        for(int i=0; i<root2d; i++) {
            int shift = (-2*j + 2*slot)%(2*slot);
            std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);
           
            plain[2*((i*root2d+j-d+slot)%slot)] = std::make_unique<HEaaN::Message>();
            *plain[2*((i*root2d+j-d+slot)%slot)] = msg;
            plain[2*((i*root2d+j-d+slot)%slot)]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    std::cout << "Debug\n";
    */

    HEaaN::Ciphertext temp = HEaaN::Ciphertext(ctxt);

    // baby
    std::vector<uint64_t> giant_idx;
    for(int i=0; i< root2d; i++) {
        int index = ((root2d*i-d)*(d-1) + slot)%slot;
        //std::cout << 2*index << '\n';
        giant_idx.push_back(2*index);
    }

    //sort(giant_idx.begin(),giant_idx.end());

    //std::cout << "debug\n";

    std::unordered_map<int, std::unique_ptr<HEaaN::Ciphertext>> giant;

    // giant_idx가 std::vector<int>라 가정
    for (int i=0; i<giant_idx.size(); i++ ) {
        //std::cout << giant_idx[i]<< '\n';
        // 1) 새로운 Ciphertext 객체를 생성해 unique_ptr에 할당
        giant[giant_idx[i]] = std::make_unique<HEaaN::Ciphertext>(eval.getContext());
        giant[giant_idx[i]]->to(HEaaN::getCurrentCudaDevice());
        // 2) leftRotate2 호출 시, 포인터가 가리키는 객체를 참조로 전달
        if(i==0) 
            eval.leftRotate(temp, giant_idx[i], *giant[giant_idx[i]]);
        else
            eval.leftRotate(*giant[giant_idx[i-1]], 2*(d-1)*root2d, *giant[giant_idx[i]]);
            //leftRotate2(eval, *giant[giant_idx[i-1]], giant_idx[i]-giant_idx[i-1], *giant[giant_idx[i]]);
    }


    

    for(int j=0; j<root2d; j++) {
        HEaaN::Ciphertext inner(ctxt);
        inner.to(HEaaN::getCurrentCudaDevice());
        bool start = true;

        for(int i=0; i<root2d; i++) {

            int flag = (i*root2d + j - d + slot) % slot;
            if (!((0 <= flag && flag < d) || (slot-d < flag))) {
                continue;
            }

            //int shift = (-2*j + 2*slot)%(2*slot);
            //std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            //std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            //HEaaN::Message msg(roted);
            //msg.to(HEaaN::getCurrentCudaDevice());

            int index = ((root2d*i-d)*(d-1) + slot) % slot;
            HEaaN::Ciphertext baby(temp);   
           
            //std::cout << "debug\n";
            eval.multWithoutRescale(*giant[2*index],*plain[2*((flag*(d-1)+slot)%slot)],baby);
           
            if (start) {
                inner = baby;
                start = false;
                
            } else {
                eval.add(inner,baby,inner);
            }
        }

        //break;
        
        eval.rescale(inner);

      
        leftRotate2(eval,inner,2*j*(d-1),inner);


        if (j==0) {
            ctxt_res = inner;
            //return;
        } else {
            eval.add(ctxt_res,inner,ctxt_res);
        }
    }
}

void phi(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt, HEaaN::Ciphertext &ctxt_res, int k, std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &plain) {
    int d= 128;
    int slot = d*d;

    HEaaN::Ciphertext temp1(eval.getContext());
    HEaaN::Ciphertext temp2(eval.getContext());

    eval.mult(ctxt,*plain[2*k],temp1);
    leftRotate2(eval,temp1,2*k,temp2);
    eval.sub(ctxt,temp1,ctxt_res);
    leftRotate2(eval,ctxt_res,2*(k-d+slot),ctxt_res);
    eval.add(ctxt_res,temp2,ctxt_res);
}

void ccmm(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt1, const HEaaN::Ciphertext &ctxt2,  HEaaN::Ciphertext &ctxt_res, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    HEaaN::Ciphertext A(ctxt1);
    HEaaN::Ciphertext B(ctxt2);

    sigma(eval,A,A,sigma_plain);
    tau(eval,B,B,tau_plain);
    

    HEaaN::Ciphertext temp(B);
    for(int k=0; k<d; k++) {
        HEaaN::Ciphertext inn1(A);

        if(k==0) {
            //std::cout << "debug\n";
            if(B.getLevel() > inn1.getLevel()) {
                eval.levelDown(B,inn1.getLevel(),B);
            }
            if(B.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1,B.getLevel(),inn1);
            }

            eval.tensor(inn1,B,ctxt_res);
            eval.rescale(ctxt_res);
        } else {
            phi(eval,inn1,inn1,k,phi_plain[k]);

            //std::cout << temp.getLevel() << " " << inn1.getLevel() << '\n';
            eval.leftRotate(temp,2*d,temp);

            if(temp.getLevel() > inn1.getLevel()) {
                eval.levelDown(temp,inn1.getLevel(),temp);
            }
            if(temp.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1,temp.getLevel(),inn1);
            }

            //std::cout << temp.getLevel() << " " << inn1.getLevel() << '\n';
            //std::cout << "debug\n";
            eval.tensor(inn1,temp,inn1);
            eval.rescale(inn1);

            //std::cout << "debug2\n";
            eval.add(inn1,ctxt_res,ctxt_res);;
        }
    }

    eval.relinearize(ctxt_res,ctxt_res);
}

void ccmm2(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt1, const HEaaN::Ciphertext &ctxt2,  HEaaN::Ciphertext &ctxt_res, int valid_row, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    HEaaN::Ciphertext outer(ctxt1);

    if(outer.getLevel() > ctxt2.getLevel()) {
        eval.levelDown(outer,ctxt2.getLevel(),outer);
    }

    

    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(outer);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        eval.add(outer,inner,outer);
        /*
        if (i == 0) {
			outer = ctxt1;
		} else {
			eval.add(outer,inner,outer);
		}*/
    }

    HEaaN::Ciphertext A(outer);
    HEaaN::Ciphertext B(ctxt2);

    sigma(eval,A,A,sigma_plain);
    tau(eval,B,B,tau_plain);
    
    HEaaN::Ciphertext temp(B);
    HEaaN::Ciphertext ext(B);

    
    for(int k=0; k<valid_row; k++) {
        HEaaN::Ciphertext inn1(A);
      
        if(k==0) {
            if(A.getLevel() > B.getLevel()) {
                eval.levelDown(A,B.getLevel(),A);
            }
            if(A.getLevel() < B.getLevel()) {
                eval.levelDown(B,A.getLevel(),B);
            }

            eval.tensor(A,B,ext);
            eval.rescale(ext);
        } else {
            phi(eval,inn1,inn1,k,phi_plain[k]);
            eval.leftRotate(temp,2*d,temp);
            
            if(temp.getLevel() > inn1.getLevel()) {
                eval.levelDown(temp,inn1.getLevel(),temp);
            }
            if(temp.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1,temp.getLevel(),inn1);
            }

            eval.tensor(inn1,temp,inn1);
            eval.rescale(inn1);

            eval.add(inn1,ext,ext);
        }
    }

    eval.relinearize(ext,ext);


    ctxt_res = ext;
    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(ctxt_res);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        eval.add(ctxt_res, inner,ctxt_res);
    }

}

void precompute_sigma_psi(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt1, std::vector<HEaaN::Ciphertext> &ctxt_res, int valid_row, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    HEaaN::Ciphertext outer(ctxt1);

    //if(outer.getLevel() > ctxt2.getLevel()) {
    //    eval.levelDown(outer,ctxt2.getLevel(),outer);
    //}

    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(outer);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        eval.add(outer,inner,outer);
        
    }

    HEaaN::Ciphertext A(outer);

    sigma(eval,A,A,sigma_plain);
    ctxt_res.push_back(A);


    //HEaaN::Ciphertext temp(B);
    //HEaaN::Ciphertext ext(B);
    for(int k=0; k<valid_row; k++) {
        HEaaN::Ciphertext inn1(A);
      
        if(k==0) {
            //eval.tensor(A,B,ext);
            //eval.rescale(ext);
        } else {
            phi(eval,inn1,inn1,k,phi_plain[k]);

            ctxt_res.push_back(inn1);
            //eval.leftRotate(temp,2*d,temp);
            
            //if(temp.getLevel() > inn1.getLevel()) {
            //    eval.levelDown(temp,inn1.getLevel(),temp);
            //}
            //if(temp.getLevel() < inn1.getLevel()) {
            //    eval.levelDown(inn1,temp.getLevel(),inn1);
            //}

            //eval.tensor(inn1,temp,inn1);
            //eval.rescale(inn1);

            //eval.add(inn1,ext,ext);
        }
    }

}

void precompute_tau_phi(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt1, std::vector<HEaaN::Ciphertext> &ctxt_res, int valid_row, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    //HEaaN::Ciphertext outer(ctxt1);

    //if(outer.getLevel() > ctxt2.getLevel()) {
    //    eval.levelDown(outer,ctxt2.getLevel(),outer);
    //}

    //for(int i=0; i<log2(d/valid_row); i++) {
    //    HEaaN::Ciphertext inner(outer);
    //    leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
    //    eval.add(outer,inner,outer);
        
    //}

    HEaaN::Ciphertext A(ctxt1);
    //std::cout << "debug\n";
    tau(eval,A,A,tau_plain);
    //std::cout << "debug\n";
    ctxt_res.push_back(A);


    //HEaaN::Ciphertext temp(B);
    //HEaaN::Ciphertext ext(B);
    HEaaN::Ciphertext inn1(A);
    for(int k=0; k<valid_row; k++) {
        
      
        if(k==0) {
            //eval.tensor(A,B,ext);
            //eval.rescale(ext);
        } else {
            //phi(eval,inn1,inn1,k,phi_plain[k]);
            eval.leftRotate(inn1,2*d,inn1);
            ctxt_res.push_back(inn1);
            
            
            //if(temp.getLevel() > inn1.getLevel()) {
            //    eval.levelDown(temp,inn1.getLevel(),temp);
            //}
            //if(temp.getLevel() < inn1.getLevel()) {
            //    eval.levelDown(inn1,temp.getLevel(),inn1);
            //}

            //eval.tensor(inn1,temp,inn1);
            //eval.rescale(inn1);

            //eval.add(inn1,ext,ext);
        }
    }
    //std::cout << "debug\n";

}


void ccmm2left(const HEaaN::HomEvaluator &eval, const std::vector<HEaaN::Ciphertext> &ctxt1, const HEaaN::Ciphertext &ctxt2,  HEaaN::Ciphertext &ctxt_res, int valid_row, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    HEaaN::Ciphertext outer(ctxt1[0]);

    //if(outer.getLevel() > ctxt2.getLevel()) {
    //    eval.levelDown(outer,ctxt2.getLevel(),outer);
    //}

    //for(int i=0; i<log2(d/valid_row); i++) {
        //HEaaN::Ciphertext inner(outer);
        //leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        //eval.add(outer,inner,outer);
        /*
        if (i == 0) {
			outer = ctxt1;
		} else {
			eval.add(outer,inner,outer);
		}*/
    //}

    HEaaN::Ciphertext A(outer);
    HEaaN::Ciphertext B(ctxt2);

    //sigma(eval,A,A,sigma_plain);
    tau(eval,B,B,tau_plain);
    
    HEaaN::Ciphertext temp(B);
    HEaaN::Ciphertext ext(B);
    for(int k=0; k<valid_row; k++) {
        HEaaN::Ciphertext inn1(ctxt1[k]);
      
        if(k==0) {
            if(temp.getLevel() > inn1.getLevel()) {
                eval.levelDown(temp,inn1.getLevel(),temp);
            }

            if(temp.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1, temp.getLevel(), inn1);
            }
            eval.tensor(inn1,temp,ext);
            eval.rescale(ext);
        } else {
            //phi(eval,inn1,inn1,k,phi_plain[k]);
            eval.leftRotate(temp,2*d,temp);
            
            if(temp.getLevel() > inn1.getLevel()) {
                eval.levelDown(temp,inn1.getLevel(),temp);
            }

            if(temp.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1, temp.getLevel(), inn1);
            }

            eval.tensor(inn1,temp,inn1);
            eval.rescale(inn1);

            eval.add(inn1,ext,ext);
        }
    }

    eval.relinearize(ext,ext);


    ctxt_res = ext;
    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(ctxt_res);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        eval.add(ctxt_res, inner,ctxt_res);
    }

}


void ccmm2right(const HEaaN::HomEvaluator &eval,  const HEaaN::Ciphertext &ctxt1,const std::vector<HEaaN::Ciphertext> &ctxt2,  HEaaN::Ciphertext &ctxt_res, int valid_row, 
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
        std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
        std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {


    int d = 128;
    HEaaN::Ciphertext outer(ctxt1);

    //if(outer.getLevel() > ctxt2[1].getLevel()) {
    //    eval.levelDown(outer,ctxt2[1].getLevel(),outer);
    //}
   

    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(outer);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        //std::cout << "debug\n";
        eval.add(outer,inner,outer);
        
        //if (i == 0) {
		//	outer = ctxt1;
		//} else {
		//	eval.add(outer,inner,outer);
		//}
    }
    //std::cout << "debug\n";

    HEaaN::Ciphertext A(outer);
    HEaaN::Ciphertext B(ctxt2[0]);

    sigma(eval,A,A,sigma_plain);

   
    //tau(eval,B,B,tau_plain);
    
    HEaaN::Ciphertext temp(B);
    HEaaN::Ciphertext ext(B);
    for(int k=0; k<valid_row; k++) {
        HEaaN::Ciphertext inn1(A);
        HEaaN::Ciphertext inn2(ctxt2[k]);
      
        if(k==0) {
            if(A.getLevel() > ctxt2[0].getLevel()) {
                eval.levelDown(A,ctxt2[0].getLevel(),A);
            }
            if(A.getLevel() < B.getLevel()) {
                eval.levelDown(B,A.getLevel(),B);
            }
            //std::cout << "debug\n";
            eval.tensor(A,B,ext);
            //std::cout << "debug\n";
            eval.rescale(ext);
            //std::cout << "debug\n";
            
        } else {
            phi(eval,inn1,inn1,k,phi_plain[k]);
            //eval.leftRotate(temp,2*d,temp);
            
            if(inn2.getLevel() > inn1.getLevel()) {
                eval.levelDown(inn2,inn1.getLevel(),inn2);
            }

            if(inn2.getLevel() < inn1.getLevel()) {
                eval.levelDown(inn1, inn2.getLevel(), inn1);
            }

            eval.tensor(inn1,inn2,inn1);
            eval.rescale(inn1);
            //std::cout << "debug\n";
            eval.add(inn1,ext,ext);
        }
    }

    eval.relinearize(ext,ext);


    ctxt_res = ext;
    for(int i=0; i<log2(d/valid_row); i++) {
        HEaaN::Ciphertext inner(ctxt_res);
        leftRotate2(eval,inner,2*(1<<i)*d*valid_row,inner);
        eval.add(ctxt_res, inner,ctxt_res);
    }

}
