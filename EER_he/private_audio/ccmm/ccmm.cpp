////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2024 CryptoLab, Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of CryptoLab, Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HEaaNTimer.hpp"
#include "HEaaN-math/HEaaN-math.hpp"
#include "examples.hpp"
#include "ccmm.hpp"
/*
In this example, we describe performing HEaaN on GPU. The following code
performs the same operations as the code in '5-bootstrap.cpp'.
*/
int main() {
    printCopyright();
    if (!HEaaN::CudaTools::isAvailable()) {
        std::cout << "Device is not available." << std::endl;
        return -1;
    }

    HEaaN::HEaaNTimer timer(true);
    // You can use other bootstrappable parameter instead of FGb.
    // See 'include/HEaaN/ParameterPreset.hpp' for more details.
    HEaaN::ParameterPreset preset = HEaaN::ParameterPreset::FGb;
    HEaaN::Context context = makeContext(preset, {0}); // set CUDA device ID
    std::cout << "Parameter : " << presetNamer(preset) << std::endl
              << std::endl;

    const auto log_slots = getLogFullSlots(context);

    HEaaN::SecretKey sk(context);
    HEaaN::KeyPack pack(context);
    HEaaN::KeyGenerator keygen(context, sk, pack);

    std::cout << "Generate encryption key ... " << std::endl;
    keygen.genEncKey();
    std::cout << "done" << std::endl << std::endl;

    std::cout << "Generate multiplication key and conjugation key ... "
              << std::endl;
    keygen.genMultKey();
    keygen.genConjKey();
    keygen.genRotKeyBundle();

    std::cout << "done" << std::endl << std::endl;

    std::cout << "Generate rotation keys used in the bootstrap process ..."
              << std::endl;
    keygen.genRotKeysForBootstrap(log_slots);
    std::cout << "done" << std::endl << std::endl;

    std::cout << "Load all public keys to GPU memory ..." << std::endl;
    timer.start("* ");
    pack.to(HEaaN::getCurrentCudaDevice());
    timer.end();
    std::cout << std::endl;

    HEaaN::Encryptor enc(context);
    HEaaN::Decryptor dec(context);
    HEaaN::EnDecoder ecd(context);
    HEaaN::HomEvaluator eval(context, pack);

    /*
    Bootstrapper constructor pre-compute the constants for bootstrapping.
    */
    std::cout << "Generate Bootstrapper (including pre-computing constants for "
                 "bootstrapping) ..."
              << std::endl;
    timer.start("* ");
    HEaaN::Bootstrapper btp(eval);
    timer.end();
    std::cout << std::endl;

    std::cout << "Load boot constants data to GPU memory ..." << std::endl;
    timer.start("* ");
    btp.loadBootConstants(log_slots, HEaaN::getCurrentCudaDevice());
    timer.end();

    // plain_matrix 
    std::unordered_map<int, std::vector<std::complex<double>>> ext_sigma = loadBinaryComplex("complex_sigma.bin");
    int d = 128; // Hard coding
    int root2d = int(ceil(sqrt(2*d)));
    int slot = d*d;

    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> sigma_plain;
    for(int j=0; j<root2d; j++) {
        for(int i=0; i<root2d; i++) {
            if (!(-d<i*root2d+j-d && i*root2d+j-d<d)) {
                continue;
            }

            int shift = (-2*j + 2*slot)%(2*slot);
            std::vector<std::complex<double>> roted = ext_sigma[2*((i*root2d+j-d+slot)%slot)];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);

            //std::cout << roted.size() << '\n';
            //std::cout << roted.size() << '\n';
            //printMessage(msg);
            //std::cout << "debug" << '\n';
            HEaaN::Plaintext ptxt = ecd.encode(msg);
            //std::cout << "debug" << '\n';
           
            sigma_plain[2*((i*root2d+j-d+slot)%slot)] = std::make_unique<HEaaN::Plaintext>(context);
            *sigma_plain[2*((i*root2d+j-d+slot)%slot)] = ptxt;
            sigma_plain[2*((i*root2d+j-d+slot)%slot)]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    // plain_matrix 
    std::unordered_map<int, std::vector<std::complex<double>>> ext_tau = loadBinaryComplex("complex_tau.bin");
    //int d = 128; // Hard coding
    int rootd = int(ceil(sqrt(d)));
    //int slot = d*d;

    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> tau_plain;
    for(int j=0; j<rootd; j++) {
        for(int i=0; i<rootd; i++) {
            
            if (!(((i*rootd+j)*d)%d == 0 && (i*rootd+j)*d < slot)) {
                continue;
            }


            int shift = (-2*j*d + 2*slot)%(2*slot);
            std::vector<std::complex<double>> roted = ext_tau[2*(i*rootd+j)*d];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);

            //std::cout << roted.size() << '\n';
            //std::cout << roted.size() << '\n';
            //printMessage(msg);
            //std::cout << "debug" << '\n';
            HEaaN::Plaintext ptxt = ecd.encode(msg);
            //std::cout << "debug" << '\n';
        
            tau_plain[2*(i*rootd+j)*d] = std::make_unique<HEaaN::Plaintext>(context);
            *tau_plain[2*(i*rootd+j)*d] = ptxt;
            tau_plain[2*(i*rootd+j)*d]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    std::unordered_map<int, std::vector<std::complex<double>>> ext_transpose = loadBinaryComplex("complex_transpose.bin");
    //int d = 128; // Hard coding
    //int rootd = int(ceil(sqrt(d)));
    //int slot = d*d;

    //int count = 0;
    //std::cout << d <<" " <<  root2d << " " << slot << '\n';
    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> transpose_plain;
    for(int j=0; j<root2d; j++) {
        for(int i=0; i<root2d; i++) {
            
            int flag = (i*root2d + j - d+slot)%slot ;
            if (!((0<=flag && flag<d ) || (slot-d<flag))) {
                continue;
            }
           
            int shift = (-2*j*(d-1) + 2*slot)%(2*slot);
     
            
            std::vector<std::complex<double>> roted = ext_transpose[2*((flag*(d-1)+slot)%slot)];
            std::rotate(roted.begin(),roted.begin()+shift,roted.end());
            HEaaN::Message msg(roted);

            HEaaN::Plaintext ptxt = ecd.encode(msg);
          
        
            transpose_plain[2*((flag*(d-1)+slot)%slot)] = std::make_unique<HEaaN::Plaintext>(context);
            *transpose_plain[2*((flag*(d-1)+slot)%slot)] = ptxt;
            transpose_plain[2*((flag*(d-1)+slot)%slot)]->to(HEaaN::getCurrentCudaDevice());
        }
    }

    std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> phi_plain(d);
    std::vector<std::unordered_map<int, std::vector<std::complex<double>>>> ext_phi = loadBinary2Complex("complex_phi.bin");

    for(int k=0; k<d; ++k) {
        std::vector<std::complex<double>> roted = ext_phi[k][2*k];
        //std::cout << "debug\n";
        std::rotate(roted.begin(),roted.begin()+2*slot-2*k,roted.end());
        //std::cout << "debug\n";
        HEaaN::Message msg(roted);
        HEaaN::Plaintext ptxt = ecd.encode(msg);

        phi_plain[k][2*k] = std::make_unique<HEaaN::Plaintext>(context);
        *phi_plain[k][2*k] = ptxt;
        phi_plain[k][2*k]->to(HEaaN::getCurrentCudaDevice());
    }
    
 
    //std::cout << "count : " << count << '\n';



    {
        using namespace std;
        cout << (1<<log_slots) << '\n';
        vector<vector<double>> mat = gen_simple_mat(128,4,4);
        vector<vector<double>> mat2 = gen_simple_mat(128,4,4);
        HEaaN::Message msg2(log_slots);
        print_mat(mat,6,6);
        print_mat(mat2,6,6);
        packed(mat,mat2,msg2);

        std::cout << std::endl << "Input message : " << std::endl;
        printMessage(msg2,false,8,8);
        std::cout << std::endl;
        //HEaaN::Message msg(log_slots);
        //fillRandomReal(msg,1<<log_slots);

     

        /*
        Use 'to' function for loading object to GPU memory.
        For all HEaaN function, if input data is on GPU then output data is also
        on GPU.
        */
        msg2.to(HEaaN::getCurrentCudaDevice());

        std::cout << "Encrypt ... ";
        HEaaN::Ciphertext ctxt(context);
        //ctxt.to(HEaaN::getCurrentCudaDevice());
        enc.encrypt(msg2, pack, ctxt);
        std::cout << "done" << std::endl << std::endl;

        std::cout << "Input ciphertext - level " << ctxt.getLevel()
                  << std::endl;


        HEaaN::Ciphertext ctxt_rst(context);
        std::vector<HEaaN::Ciphertext> precom;
        //precompute_sigma_psi(eval,ctxt,precom,32,sigma_plain,tau_plain,phi_plain);
        precompute_tau_phi(eval,ctxt,precom,32,sigma_plain,tau_plain,phi_plain);
        timer.start("* ");
        //leftRotate2(eval,ctxt,5,ctxt_rst);
        //sigma(eval,ctxt,ctxt_rst,sigma_plain);
        //tau(eval,ctxt,ctxt_rst,tau_plain);
        //transpose(eval,ctxt,ctxt_rst,transpose_plain);
        //for(int i=0; i<d; i++){
        //    phi(eval,ctxt,ctxt_rst,i,phi_plain[i]);
        //}
        std::cout << ctxt.getLevel()<< '\n';
        //ccmm(eval,ctxt,ctxt,ctxt_rst,sigma_plain,tau_plain,phi_plain);
        //ccmm2(eval,ctxt,ctxt,ctxt_rst,32,sigma_plain,tau_plain,phi_plain);
        //ccmm2left(eval,precom,ctxt,ctxt_rst,32,sigma_plain,tau_plain,phi_plain);
        ccmm2right(eval,ctxt,precom,ctxt_rst,32,sigma_plain,tau_plain,phi_plain);
        std::cout << ctxt_rst.getLevel()<< '\n';
        timer.end();
        /*
        Decrypt does not work well
        if ciphertext and key data are on different device memory.
        */
        //std::cout << "debug\n";
        sk.to(HEaaN::getCurrentCudaDevice());
        ctxt_rst.to(HEaaN::getCurrentCudaDevice());
        //std::cout << "debug\n";

        HEaaN::Message dmsg, dmsg_boot;
        std::cout << "Decrypt ... ";
        dec.decrypt(ctxt_rst, sk, dmsg);
        //dec.decrypt(ctxt_boot, sk, dmsg_boot);
        std::cout << "done" << std::endl;

        /*
        Decrypted message data should be on CPU memory
        because below code is not working for GPU
        */
        dmsg.to(HEaaN::getDefaultDevice());
        dmsg_boot.to(HEaaN::getDefaultDevice());

            
        unpacked(mat,mat2,dmsg);
        print_mat(mat,6,6);
        print_mat(mat2,6,6);

        std::cout.precision(10);
        std::cout << std::endl << "Decrypted message : " << std::endl;
        printMessage(dmsg,8,8);
    }

    return 0;
}
