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
#include "activation.hpp"

void approximate_derivated_ReLU(const HEaaN::HomEvaluator &eval,
                                const HEaaN::Bootstrapper &btp,
                                const HEaaN::Ciphertext &ctxt,
                                HEaaN::Ciphertext &ctxt_res,
                                const HEaaN::Real multiplier);

void approximate_ReLU(const HEaaN::HomEvaluator &eval,
                      const HEaaN::Bootstrapper &btp,
                      const HEaaN::Ciphertext &ctxt,
                      HEaaN::Ciphertext &ctxt_res,
                      const HEaaN::Real multiplier);

void approximate_Sigmoid(const HEaaN::HomEvaluator &eval,
                        const HEaaN::Bootstrapper &btp,
                        const HEaaN::Ciphertext &ctxt,
                        HEaaN::Ciphertext &ctxt_res,
                        const HEaaN::Real multiplier);
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

    {
        HEaaN::Message msg(log_slots);
        fillRandomReal(msg,1<<log_slots);

        std::cout << std::endl << "Input message : " << std::endl;
        printMessage(msg);
        std::cout << std::endl;

        /*
        Use 'to' function for loading object to GPU memory.
        For all HEaaN function, if input data is on GPU then output data is also
        on GPU.
        */
        msg.to(HEaaN::getCurrentCudaDevice());

        std::cout << "Encrypt ... ";
        HEaaN::Ciphertext ctxt(context);
        enc.encrypt(msg, pack, ctxt);
        std::cout << "done" << std::endl << std::endl;

        std::cout << "Input ciphertext - level " << ctxt.getLevel()
                  << std::endl;


        HEaaN::Ciphertext ctxt_rst(context);

        timer.start("* ");
        //approximate_derivated_ReLU(eval,btp,ctxt,ctxt_rst,1.0);
        //approximate_ReLU(eval,btp,ctxt,ctxt_rst,1.0);
        approximate_Sigmoid(eval,btp,ctxt,ctxt_rst,1.0);
        timer.end();
        /*
        Decrypt does not work well
        if ciphertext and key data are on different device memory.
        */
        sk.to(HEaaN::getCurrentCudaDevice());

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

        std::cout.precision(10);
        std::cout << std::endl << "Decrypted message : " << std::endl;
        printMessage(dmsg);
    }

    return 0;
}


