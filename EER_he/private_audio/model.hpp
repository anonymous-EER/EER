#pragma once

//#include "HEaaNTimer.hpp"
#include "HEaaN-math/HEaaN-math.hpp"
//#include "examples.hpp"
#include "ccmm/ccmm.hpp"
#include "activation/activation.hpp"

using namespace HEaaN;


class Model {
public :
    // training
    int n_sample;
    int epoch;
    double lr;
    int input_dim1, input_dim2;
    int hidden_dim;
    int log_slot = 15;

    std::unique_ptr<HEaaN::Message> add_bias ;
    std::unique_ptr<HEaaN::Message> add_bias2 ;

    Model(int _n, int _epoch, double _lr, int _dim1, int _dim2, int _hdim) {
        n_sample = _n;
        epoch = _epoch;
        lr = _lr;
        input_dim1 = _dim1;
        input_dim2 = _dim2;
        hidden_dim = _hdim;
        add_bias = std::make_unique<HEaaN::Message>(log_slot,std::complex<double>(0.0));
        add_bias2 = std::make_unique<HEaaN::Message>(log_slot,std::complex<double>(0.0));
    };


    void create_masking(const HomEvaluator &eval) {
      
        std::vector<std::vector<double>> mat(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i < this->n_sample && j==this->input_dim1)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }

        std::vector<std::vector<double>> mat2(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i < this->n_sample && j==this->input_dim2)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,mask);
        //printMessage(mask);
        *add_bias = mask;
        add_bias->to(HEaaN::getCurrentCudaDevice());

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j< this->n_sample && i ==this->hidden_dim)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->n_sample && i==this->hidden_dim)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        packed(mat,mat2,mask);
        //printMessage(mask);
        *add_bias2 = mask;
        add_bias2->to(HEaaN::getCurrentCudaDevice());

        
    }

   

    void inference(const HomEvaluator &eval,const Bootstrapper &btp,const std::vector<Ciphertext> &precom_w1,const std::vector<Ciphertext> &precom_w2, Ciphertext& x, Ciphertext & y,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
                    std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {
        

        Ciphertext z1(eval.getContext());
        //std::cout << "debug\n";

        if(x.getLevel() < 3+btp.getMinLevelForBootstrap()) {
            std::cout << "debug\n";
            btp.bootstrapExtended(x,x,false);
        }

 

        ccmm2left(eval,precom_w1,x,z1,32,sigma_plain,tau_plain,phi_plain);
       
        Ciphertext a1(eval.getContext());
        approximate_ReLU(eval,btp,z1,a1,1.0);
       

        eval.add(a1,*add_bias2,a1);
                    
        Ciphertext z2(eval.getContext());


        if(a1.getLevel() < 3+btp.getMinLevelForBootstrap()) {
            std::cout << " ?\n";
           
        }
        btp.bootstrapExtended(a1,a1,false);

        ccmm2left(eval,precom_w2,a1,z2,1,sigma_plain,tau_plain,phi_plain);

        //y=Ciphertext(z2);
        //return;

        approximate_Sigmoid(eval,btp,z2,y,1.0);
    }

    void forward(const HomEvaluator &eval,const Bootstrapper &btp,const Ciphertext &w1,const Ciphertext &w2, std::vector<Ciphertext>& precom_x, 
                    Ciphertext & y, Ciphertext &sign_z1, Ciphertext &a1_temp, HEaaN::Message &sigmoid_mask,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
                    std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {
        

        Ciphertext z1(eval.getContext());
        //std::cout << "debug\n";

        //if(x.getLevel() < 3+btp.getMinLevelForBootstrap()) {
        //    std::cout << "debug\n";
        //    btp.bootstrap(x,x,false);
        //}

 

        ccmm2right(eval,w1,precom_x,z1,32,sigma_plain,tau_plain,phi_plain);
       
        Ciphertext a1(eval.getContext());

        if(z1.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
            //std::cout << "btp\n";
            btp.bootstrapExtended(z1,z1,false);
        }

        // scaling 
        HEaaN::Ciphertext temp = HEaaN::Ciphertext(z1);

        eval.mult(temp,0.01,temp);

        approximate_derivated_ReLU(eval,btp,temp,a1,1.0);

        //save
        sign_z1 = Ciphertext(a1);

        eval.mult(z1,a1,a1);
        btp.bootstrapExtended(a1,a1,false);

        eval.add(a1,*add_bias2,a1);

        //save
        a1_temp = Ciphertext(a1);

        Ciphertext z2(eval.getContext());


        if(a1.getLevel() < 3+btp.getMinLevelForBootstrap()) {
            //std::cout << " ?\n";
            btp.bootstrapExtended(a1,a1,false);
        }

        
        ccmm2(eval,w2,a1,z2,1,sigma_plain,tau_plain,phi_plain);
        approximate_Sigmoid(eval,btp,z2,y,1.0);
        //y = Ciphertext(a1);
        btp.bootstrap(y,y,false);

        eval.mult(y,sigmoid_mask,y);
    }

    void backward(const Decryptor &dec, SecretKey &sk, HomEvaluator &eval,const Bootstrapper &btp, Ciphertext &w1, Ciphertext &w2, std::vector<Ciphertext>& precom_x, Ciphertext & y,
                    Ciphertext &label,Ciphertext &onehot, Ciphertext &sign_z1, Ciphertext &a1,
                    Message &del_bias_mask,Message &sign_mask,
                    Ciphertext &dLdw2, Ciphertext &dLdw1,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &transpose_plain,
                    std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {

        Ciphertext dLdz(eval.getContext());
        //std::cout << "debug1\n";
        //btp.bootstrap(y,y);
        //std::cout << "debug2\n";
          //std::cout << "debug\n";
        eval.sub(y,label,dLdz);
          //std::cout << "debug\n";

          


        eval.mult(dLdz,onehot,dLdz);

        
        std::vector<std::vector<double>> mat(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->n_sample && i<this->hidden_dim+1)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }

        std::vector<std::vector<double>> mat2(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->n_sample && i<this->hidden_dim+1)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,mask);
        //printMessage(mask);
        //*add_bias = mask;
        mask.to(HEaaN::getCurrentCudaDevice());
        //std::cout << "debug\n";
        eval.mult(a1,mask,a1);
        //std::cout << "debug\n";
        //Ciphertext dLdw2(eval.getContext());
        Ciphertext a1_temp(a1);
        //std::cout << "debug\n";
        transpose(eval,a1,a1_temp,transpose_plain);
        
        //std::cout << dLdz.getLevel() << '\n';

        //btp.bootstrap(a1_temp,a1_temp,false);
        //std::cout << a1_temp.getLevel() << '\n';
        ccmm2(eval,dLdz,a1_temp,dLdw2,1,sigma_plain,tau_plain,phi_plain);
        //std::cout << "debug\n";
         
        eval.mult(dLdw2,1.0/this->n_sample,dLdw2);

       
        //std::cout << "debug1\n";
        //dLdw2 = Ciphertext(dLdz);
        //std::cout << "debug2\n";
        //return;

        Ciphertext trans_w2(w2);
        //std::cout << "debug\n";
        eval.mult(trans_w2,del_bias_mask,trans_w2);
        //std::cout << "debug\n";
        transpose(eval,w2,trans_w2,transpose_plain);
        //std::cout << "debug\n";
        Ciphertext temp(eval.getContext());

        ccmm(eval,trans_w2,dLdz,temp,sigma_plain,tau_plain,phi_plain);

        //std::cout << "debug\n";
        eval.mult(temp,sign_z1,temp);
        //std::cout << "debug\n";
        eval.mult(temp,sign_mask,temp);
        //std::cout << "debug\n";
        
        
        btp.bootstrapExtended(temp,temp,false);
        ccmm2right(eval,temp,precom_x,dLdw1,32,sigma_plain,tau_plain,phi_plain);
        eval.mult(dLdw1,1.0/this->n_sample,dLdw1);

        //std::cout << "debug\n";

        


    }

    void cotrain(const Decryptor&dec, SecretKey &sk, HomEvaluator &eval,const Bootstrapper &btp, Ciphertext &w1, Ciphertext &w2, Ciphertext& x, Ciphertext & y, int round,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &sigma_plain,
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &tau_plain, 
                    std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>> &transpose_plain,
                    std::vector<std::unordered_map<int,std::unique_ptr<HEaaN::Plaintext>>> &phi_plain) {

        std::vector<std::vector<double>> mat(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i==0 && j<this->n_sample)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }

        std::vector<std::vector<double>> mat2(128,std::vector<double>(128,0));

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i==0 && j<this->n_sample)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message sigmoid_mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,sigmoid_mask);
        sigmoid_mask.to(HEaaN::getCurrentCudaDevice());



        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->input_dim1+1 && i < this->hidden_dim)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->input_dim2+1 && i < this->hidden_dim)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message weight_mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,weight_mask);
        weight_mask.to(HEaaN::getCurrentCudaDevice());

        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->hidden_dim+1 && i < 1)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j < this->hidden_dim+1 && i < 1)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message weight_mask2(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,weight_mask2);
        weight_mask2.to(HEaaN::getCurrentCudaDevice());


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j< this->n_sample && i < this->hidden_dim)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (j< this->n_sample && i < this->hidden_dim)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message sign_mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,sign_mask);
        sign_mask.to(HEaaN::getCurrentCudaDevice());


         for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i < 128 && j < this->hidden_dim)
                    mat[i][j] = 1.0;
                else 
                    mat[i][j] = 0.0;
            }
        }


        for(int i=0; i<128; i++) {
            for(int j=0; j<128; j++) {
                if (i < 128 && j < this->hidden_dim)
                    mat2[i][j] = 1.0;
                else 
                    mat2[i][j] = 0.0;
            }
        }

       
        HEaaN::Message del_bias_mask(this->log_slot,std::complex<double>(0.0));
        packed(mat,mat2,del_bias_mask);
        del_bias_mask.to(HEaaN::getCurrentCudaDevice());


        std::vector<std::complex<double>> vec((1<<log_slot),0.0);

        for(int i=0; i<vec.size()/2; i++) {
            vec[2*i] = 1.0;
            vec[2*i+1] = 0.0;
        }

        HEaaN::Message mask(vec);
        mask.to(HEaaN::getCurrentCudaDevice());


        // cococococococo
        //eval.add(x,*add_bias,x);
        //transpose(eval,x,x,transpose_plain);

        
        std::vector<Ciphertext> precom_x;
        precompute_tau_phi(eval,x,precom_x,32,sigma_plain,tau_plain,phi_plain);
        
        Ciphertext x_T(x);
        transpose(eval,x,x_T,transpose_plain);
        std::vector<Ciphertext> precom_xT;
        precompute_tau_phi(eval,x_T,precom_xT,32,sigma_plain,tau_plain,phi_plain);

      
        Ciphertext sign_z1(eval.getContext()), a1(eval.getContext());

        for(int rd =0 ;  rd<round; rd++) {
            forward(eval,btp,w1,w2,precom_x,y,sign_z1,a1,sigmoid_mask,sigma_plain,tau_plain,phi_plain);
            //std::cout << "Debug\n";

            
            //Message dmsg;
            //sk.to(getCurrentCudaDevice());
            //dec.decrypt(y,sk,dmsg);
            //dmsg.to(HEaaN::getDefaultDevice());

                


            //std::cout << "a1 ... \n";
            //printMessage(dmsg,false,8);

            HEaaN::Ciphertext label(eval.getContext());
            HEaaN::Ciphertext label_cp(eval.getContext());

            eval.sub(y,0.5,label);
            label_cp = Ciphertext(label);
            eval.leftRotate(label_cp,1,label_cp);

            eval.mult(label,mask,label);
            eval.mult(label_cp,mask,label_cp);

            eval.leftRotate(label,-1,label);

            eval.add(label,label_cp,label);

            approximate_derivated_ReLU(eval,btp,label,label,1.0);
            eval.mult(label,sigmoid_mask,label);
            btp.bootstrap(label,label,false);

            HEaaN::Ciphertext tr(eval.getContext());
            HEaaN::Ciphertext tr_cp(eval.getContext());

            eval.sub(y,0.5,tr);
            eval.square(tr,tr);
            eval.sub(tr,0.09,tr);

            tr_cp = Ciphertext(tr);
            eval.leftRotate(tr_cp,1,tr_cp);

            eval.mult(tr,mask,tr);
            eval.mult(tr_cp,mask,tr_cp);

            eval.leftRotate(tr,-1,tr);

            eval.add(tr,tr_cp,tr);

            HEaaN::Ciphertext onehot(eval.getContext());
            approximate_derivated_ReLU(eval,btp,tr,onehot,1.0);
            eval.mult(onehot,sigmoid_mask,onehot);
            btp.bootstrap(onehot,onehot,false);


            
            /*Message dmsg;
            sk.to(getCurrentCudaDevice());
            dec.decrypt(label,sk,dmsg);
            dmsg.to(HEaaN::getDefaultDevice());

            std::cout << "w1 ... \n";
            printMessage(dmsg,false,8);

            

            sk.to(getCurrentCudaDevice());
            dec.decrypt(onehot,sk,dmsg);
            dmsg.to(HEaaN::getDefaultDevice());

            std::cout << "w2 ... \n";
            printMessage(dmsg,false,8);*/


            for(int i=0; i<this->epoch; i++) {
                
                HEaaN::HEaaNTimer timer(true);
                timer.start("* epoch");
                Ciphertext sign_z1(eval.getContext()), a1(eval.getContext());
                Ciphertext dLdw1(eval.getContext()), dLdw2(eval.getContext());
                forward(eval,btp,w1,w2,precom_x,y,sign_z1,a1,sigmoid_mask,sigma_plain,tau_plain,phi_plain);

               
                
                backward(dec,sk,eval,btp,w1,w2,precom_xT,y,label,onehot,sign_z1,a1,del_bias_mask,sign_mask,dLdw2,dLdw1,sigma_plain,tau_plain,transpose_plain,phi_plain);
                
               


                eval.mult(dLdw1,this->lr,dLdw1);
                eval.sub(w1,dLdw1,w1);
                eval.mult(dLdw2,this->lr,dLdw2);
                eval.sub(w2,dLdw2,w2);


                eval.mult(w1,weight_mask,w1);
                eval.mult(w2,weight_mask2,w2);

                


                //std::cout << w1.getLevel() << '\n';
                //std::cout << w2.getLevel() << '\n';
                btp.bootstrapExtended(w1,w1,false);
                btp.bootstrapExtended(w2,w2,false);
                timer.end();
            }
        }

    }

};