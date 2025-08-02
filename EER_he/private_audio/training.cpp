

#include "HEaaN-math/HEaaN-math.hpp"
#include "ccmm/ccmm.hpp"
#include "activation/activation.hpp"
#include "json.hpp" // 직접 다운로드한 헤더
#include "model.hpp"

using json = nlohmann::json;

std::vector<std::vector<std::vector<double>>> split_matrix(
    const std::vector<std::vector<double>>& matrix,
    size_t chunk_rows
) {
    std::vector<std::vector<std::vector<double>>> chunks;
    size_t total_rows = matrix.size();
    size_t cols = matrix[0].size();

    for (size_t i = 0; i < total_rows; i += chunk_rows) {
        size_t end_row = std::min(i + chunk_rows, total_rows);
        std::vector<std::vector<double>> chunk;

        for (size_t j = i; j < end_row; ++j) {
            chunk.push_back(matrix[j]);
        }

        chunks.push_back(chunk);
    }

    return chunks;
}


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

    std::vector<double> accs;

    for(int pnum = 12; pnum < 13; pnum++) {
        std::ifstream file("model_datas_final/Given/2/"+std::to_string(pnum)+"_100_0.json");
        //std::ifstream file("cotrain/model_saves/Wesad/" + std::to_string(pnum) + "_100.json");
        if (!file.is_open()) {
            std::cerr << "JSON 파일을 열 수 없습니다.\n";
            return 1;
        }

        json j;
        file >> j;

        std::vector<std::vector<std::vector<double>>> x_preds;
        std::vector<std::vector<std::vector<double>>> x_tests;
        //std::vector<double> y_test;
        
        for (const auto& sample : j["x_preds"]) {
            std::vector<std::vector<double>> sample_vec;
            for (const auto& row : sample) {
                sample_vec.push_back(row.get<std::vector<double>>());
            }
            x_preds.push_back(sample_vec);
        }

       // std::cout << "x_preds 크기: " << x_preds.size() << " x "
         //      << x_preds[0].size() << " x "
          //      << x_preds[0][0].size() << " & " << x_preds[1][0].size() <<  std::endl;

        for (const auto& sample : j["x_tests"]) {
            std::vector<std::vector<double>> sample_vec;
            for (const auto& row : sample) {
                sample_vec.push_back(row.get<std::vector<double>>());
            }
            x_tests.push_back(sample_vec);
        }
                
        //std::cout << "x_tests 크기: " << x_tests.size() << " x "
         //           << x_tests[0].size() << " x "
         //           << x_tests[0][0].size() << std::endl;

        std::vector<std::vector<double>> y_test;

        if (j.contains("y_test")) {
            y_test = j["y_test"].get<std::vector<std::vector<double>>>();
        } else {
            std::cerr << "\"y_test\" 항목이 JSON에 없습니다.\n";
        }

       // std::cout << "y_test 크기: " << y_test.size() << " x "
        //            << y_test[0].size()  << std::endl;


        
        file.close();
        
        std::ifstream file2("model_saves_final/Given/2/"+std::to_string(pnum)+"_100_0_0.json");
        //std::ifstream file2("cotrain/model_datas/Wesad/" + std::to_string(pnum) + "_0.json");
        if (!file2.is_open()) {
            std::cerr << "파일 열기 실패\n";
            return 1;
        }

        file2 >> j;

        // w1, w2 변수 선언
        std::vector<std::vector<double>> w1_0;
        std::vector<std::vector<double>> w2_0;

        try {
            if (j.contains("w1")) {
                w1_0 = j["w1"].get<std::vector<std::vector<double>>>();
                //std::cout << "w1 크기: " << w1_0.size() << " x " 
                //        << (w1_0.empty() ? 0 : w1_0[0].size()) << std::endl;
            } else {
                std::cerr << "[경고] 'w1' 키가 JSON에 없습니다.\n";
            }

            if (j.contains("w2")) {
                w2_0 = j["w2"].get<std::vector<std::vector<double>>>();
                //std::cout << "w2 크기: " << w2_0.size() << " x " 
                //        << (w2_0.empty() ? 0 : w2_0[0].size()) << std::endl;
            } else {
                std::cerr << "[경고] 'w2' 키가 JSON에 없습니다.\n";
            }

        } catch (const json::exception& e) {
            std::cerr << "[에러] JSON 파싱 실패: " << e.what() << std::endl;
            return 1;
        }

        // 예시 출력
        if (!w1_0.empty() && !w1_0[0].empty()) {
            //std::cout << "w1[0][0] = " << w1_0[0][0] << std::endl;
        }
        if (!w2_0.empty() && !w2_0[0].empty()) {
            //std::cout << "w2[0][0] = " << w2_0[0][0] << std::endl;
        }   

        file2.close();

        std::ifstream file3("model_saves_final/Given/2/"+std::to_string(pnum)+"_100_1_0.json");
        //std::ifstream file3("cotrain/model_datas/Wesad/" + std::to_string(pnum) + "_1.json");
        if (!file3.is_open()) {
            std::cerr << "파일 열기 실패\n";
            return 1;
        }

        file3 >> j;

        // w1, w2 변수 선언
        std::vector<std::vector<double>> w1_1;
        std::vector<std::vector<double>> w2_1;

        try {
            if (j.contains("w1")) {
                w1_1 = j["w1"].get<std::vector<std::vector<double>>>();
                //std::cout << "w1 크기: " << w1_1.size() << " x " 
                 //       << (w1_1.empty() ? 0 : w1_1[0].size()) << std::endl;
            } else {
                std::cerr << "[경고] 'w1' 키가 JSON에 없습니다.\n";
            }

            if (j.contains("w2")) {
                w2_1 = j["w2"].get<std::vector<std::vector<double>>>();
                //std::cout << "w2 크기: " << w2_1.size() << " x " 
                 //       << (w2_1.empty() ? 0 : w2_1[0].size()) << std::endl;
            } else {
                std::cerr << "[경고] 'w2' 키가 JSON에 없습니다.\n";
            }

        } catch (const json::exception& e) {
            std::cerr << "[에러] JSON 파싱 실패: " << e.what() << std::endl;
            return 1;
        }

        // 예시 출력
        if (!w1_1.empty() && !w1_1[0].empty()) {
            //std::cout << "w1[0][0] = " << w1_1[0][0] << std::endl;
        }
        if (!w2_1.empty() && !w2_1[0].empty()) {
            //std::cout << "w2[0][0] = " << w2_1[0][0] << std::endl;
        }

        file3.close();


        Model model(100, 15,0.001,84,127,32);
        //Model model(100, 10,0.05,61,82,32);
        model.create_masking(eval);

        HEaaN::Ciphertext ctxt_w1(context);
        HEaaN::Message packed_w1(log_slots,std::complex<double>(0.0));
        packed(w1_0,w1_1,packed_w1);
        //printMessage(packed_w1);
        packed_w1.to(HEaaN::getCurrentCudaDevice());
        enc.encrypt(packed_w1,pack,ctxt_w1);
        ctxt_w1.to(HEaaN::getCurrentCudaDevice());
        transpose(eval,ctxt_w1,ctxt_w1,transpose_plain);



        HEaaN::Message packed_w2(log_slots, std::complex<double>(0.0));
        packed(w2_0,w2_1,packed_w2);
        //printMessage(packed_w2);
        packed_w2.to(HEaaN::getCurrentCudaDevice());
        //pack.to(HEaaN::getCurrentCudaDevice());

        HEaaN::Ciphertext ctxt_w2(context);
    
        enc.encrypt(packed_w2,pack,ctxt_w2);
        ctxt_w2.to(HEaaN::getCurrentCudaDevice());
        transpose(eval,ctxt_w2,ctxt_w2,transpose_plain);
        
       


        // 분할 수행 (120 x 84 단위로)
        auto chunks = split_matrix(x_tests[0], 120);

        // 결과 출력
        //std::cout << "총 조각 개수: " << chunks.size() << std::endl;
        for (size_t i = 0; i < chunks.size(); ++i) {
            //std::cout << "Chunk " << i << " 크기: " 
             //       << chunks[i].size() << " x " 
              //      << (chunks[i].empty() ? 0 : chunks[i][0].size()) << std::endl;
        }

        auto chunks2 = split_matrix(x_tests[1], 120);

        // 결과 출력
        //std::cout << "총 조각 개수: " << chunks2.size() << std::endl;
        for (size_t i = 0; i < chunks2.size(); ++i) {
            //std::cout << "Chunk2 " << i << " 크기: " 
            //        << chunks2[i].size() << " x " 
            //        << (chunks2[i].empty() ? 0 : chunks2[i][0].size()) << std::endl;
        }


        using namespace std; 
        vector<vector<double>> pred1;
        vector<vector<double>> pred2;

        for(int i=0; i<1; i++) {
            HEaaN::Ciphertext ctxt(context);
            HEaaN::Message packed_x(log_slots,std::complex<double>(0.0,0.0));
            //print_mat(chunks[i],2,10);
            //print_mat(chunks2[i],2,10);

            packed(x_preds[0],x_preds[1],packed_x);
            //printMessage(packed_x);
            packed_x.to(HEaaN::getCurrentCudaDevice());
            //pack.to(HEaaN::getCurrentCudaDevice());
            enc.encrypt(packed_x,pack,ctxt);
        
            eval.add(ctxt,*model.add_bias,ctxt);
            transpose(eval,ctxt,ctxt,transpose_plain);


            //sk.to(HEaaN::getCurrentCudaDevice());
            //std::cout << "debug\n";
            //dec.decrypt(ctxt,sk,packed_x);
            //packed_x.to(HEaaN::getDefaultDevice());
            //printMessage(packed_x,false,8);
            //unpacked(chunks[i],chunks2[i],packed_x);
            //print_mat(chunks[i],2,10);
            //print_mat(chunks2[i],2,10);

        
            HEaaN::Ciphertext ctxt_res(context);
            HEaaN::Message packed_y;
            HEaaN::Message packed_y2;

            timer.start("Inference ");
            model.cotrain(dec,sk,eval,btp,ctxt_w1,ctxt_w2,ctxt,ctxt_res,25,sigma_plain,tau_plain,transpose_plain,phi_plain);
            timer.end();
            
         
        
            //std::cout << "debug\n";
            sk.to(HEaaN::getCurrentCudaDevice());
            //std::cout << "debug\n";
            dec.decrypt(ctxt_w1,sk,packed_y);
            //std::cout << "debug\n";
            packed_y.to(HEaaN::getDefaultDevice());
            //printMessage(packed_y);

            //std::cout << "debug\n";
            sk.to(HEaaN::getCurrentCudaDevice());
            //std::cout << "debug\n";
            dec.decrypt(ctxt_w2,sk,packed_y2);
            //std::cout << "debug\n";
            packed_y2.to(HEaaN::getDefaultDevice());
            //printMessage(packed_y);

            std::vector<std::vector<double>> y1(128,std::vector<double>(128,0.0));
            std::vector<std::vector<double>> y2(128,std::vector<double>(128,0.0));
            unpacked(y1,y2,packed_y);

            std::vector<std::vector<double>> y3(128,std::vector<double>(128,0.0));
            std::vector<std::vector<double>> y4(128,std::vector<double>(128,0.0));
            unpacked(y3,y4,packed_y2);
            
            //print_mat(y1,1,32);
            //print_mat(y2,1,32);
            json j;

            j["w1"] = y1;
            j["w2"] = y3;

            // 파일로 저장
            //std::ofstream fifile("trained_saves/Given/"+std::to_string(pnum)+"_100_0_0.json");
            std::ofstream fifile("trained_saves/Wesad/"+std::to_string(pnum)+"_100_0.json");
            if (!fifile.is_open()) {
                std::cerr << "파일 저장 실패" << std::endl;
                return 1;
            }

            // pretty print로 저장 (공백/들여쓰기 있음)
            fifile << j.dump(2);  // 인자 2는 들여쓰기 단위
            fifile.close();

            std::cout << "JSON 저장 완료: output.json" << std::endl;

            json jj;
            jj["w1"] = y2;
            jj["w2"] = y4;

            // 파일로 저장
            //std::ofstream fifile2("trained_saves/Given/"+std::to_string(pnum)+"_100_1_0.json");
            std::ofstream fifile2("trained_saves/Wesad/"+std::to_string(pnum)+"_100_1.json");
            if (!fifile2.is_open()) {
                std::cerr << "파일 저장 실패" << std::endl;
                return 1;
            }

            // pretty print로 저장 (공백/들여쓰기 있음)
            fifile2 << jj.dump(2);  // 인자 2는 들여쓰기 단위
            fifile2.close();

            std::cout << "JSON 저장 완료: output.json" << std::endl;

            //cout << "chunk size : " << chunks[i].size() << '\n';
        }

        //double acc = 0.0;
        //for(int i=0; i<y_test.size(); i++) {
         //   if(round((pred1[i][0]+pred2[i][0])/2.0) == y_test[i][0]) {
           //     acc+=1.0;
            //}
        //}

        std::cout << "pnum" << std::to_string(pnum)<< "\'s train complete\n";
        //accs.push_back(acc/y_test.size());
    }   

    //double mean = 0.0;
    //for(int i=0; i<accs.size(); i++) {
    //    mean+=accs[i];
    //}

    //std::cout << "Leave out acc : " << mean/(double)accs.size() << '\n';

    return 0;
}