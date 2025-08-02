#pragma once

#include "HEaaN/HEaaN.hpp"
#include "HEaaN-math/HEaaN-math.hpp"

void approximate_derivated_ReLU(const HEaaN::HomEvaluator &eval,
    const HEaaN::Bootstrapper &btp,
    const HEaaN::Ciphertext &ctxt,
    HEaaN::Ciphertext &ctxt_res,
    const HEaaN::Real multiplier=1.0){


std::vector<HEaaN::Real> coeffs1 = {
0.32605531770751276e-27, 0.63713355013748307, 0.96638008325953307e-27,
-0.21376474517603805, 0.14078053953548597e-26, 0.12997879454442307,
0.32384267959866018e-27, -0.0947909612869684, -0.44470882359357376e-27,
0.075917936322122959, 0.18963482898054708e-26, -0.064613710325904794,
-0.5775075296998782e-27, 0.057593917557167362, -0.31866109262040246e-26,
-0.51453436690426679
};

std::vector<HEaaN::Real> coeffs2 = {
-0.15481815164177629e-26, -0.72535634150240158, -0.30583712913994119e-26,
-1.700553029987224, -0.30226446764082583e-26, -1.253961804370583,
-0.28921966895028109e-26, -1.4441209142141685, -0.27195500541793715e-26,
-1.1666282999498775, -0.2488197447639343e-26, -1.2182643955976286,
-0.23818895236993765e-26, -0.97769563533053216, -0.21384234078281345e-26,
-0.96714980183541921, -0.18374655471146771e-26, -0.74621079317205274,
-0.15015506031038111e-26, -0.7093785137759101, -0.13423719649247374e-26,
-0.51064903419936536, -0.14860402249885922e-26, -0.47114290730816244,
-0.12948243393594419e-26, -0.30068435373287239, -0.60205084799546203e-27,
-0.4999749972462898};

std::vector<HEaaN::Real> coeffs3 = { 0.23057814159951059e-30, 1.3126801750100274, 0.45369022293103483e-30,
0.44326455937581305, 0.41987062770006384e-30, 0.7098004895996092,
0.37589409245773775e-30, 0.44185207611914302, 0.31545398476410848e-30,
0.48439982222240077, 0.25567849958108216e-30, 0.30464925217763116,
0.19270058188369575e-30, 0.28329208430576176, 0.13985081837410144e-30,
0.16401146707396689, 0.92947406067029118e-31, 0.13410380654010332,
0.59187477308104792e-31, 0.067149530888032436, 0.33506000473740497e-31,
0.048500104899959402, 0.1802868885516293e-31, 0.019142530476768867,
0.8081207422310361e-32, 0.012086934479667787, 0.33901718403134585e-32,
0.0029636098860223543, 0.96597996471713936e-33, 0.0015994759535727784 };

const HEaaN::Math::ChebyshevCoefficients cheby_coeffs1(coeffs1, 16,
HEaaN::Math::PolynomialBasis::basic);
const HEaaN::Math::ChebyshevCoefficients cheby_coeffs2(coeffs2, 32,
HEaaN::Math::PolynomialBasis::basic);
const HEaaN::Math::ChebyshevCoefficients cheby_coeffs3(coeffs3, 32,
HEaaN::Math::PolynomialBasis::basic);


//std::cout << "Level : " << ctxt_res.getLevel() << '\n';
//std::cout << btp.getMinLevelForBootstrap() << '\n';

HEaaN::Ciphertext temp = HEaaN::Ciphertext(ctxt);
if (ctxt.getLevel() < cheby_coeffs1.level_cost + btp.getMinLevelForBootstrap()) {
btp.bootstrap(ctxt,temp,false);
}

ctxt_res = HEaaN::Math::evaluateChebyshevExpansion(eval,btp,temp,cheby_coeffs1,1.0);

if (ctxt_res.getLevel() < cheby_coeffs2.level_cost + btp.getMinLevelForBootstrap()) {
btp.bootstrap(ctxt_res,ctxt_res,false);
}

ctxt_res = HEaaN::Math::evaluateChebyshevExpansion(eval,btp,ctxt_res,cheby_coeffs2,1.0);

if (ctxt_res.getLevel() < cheby_coeffs3.level_cost + btp.getMinLevelForBootstrap()) {
btp.bootstrap(ctxt_res,ctxt_res,false);
}

ctxt_res = HEaaN::Math::evaluateChebyshevExpansion(eval,btp,ctxt_res,cheby_coeffs3,1.0);


eval.add(ctxt_res, 0.5, ctxt_res);

}

void approximate_ReLU(const HEaaN::HomEvaluator &eval,
const HEaaN::Bootstrapper &btp,
HEaaN::Ciphertext &ctxt,
HEaaN::Ciphertext &ctxt_res,
const HEaaN::Real multiplier) {

if(ctxt.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
btp.bootstrap(ctxt,ctxt,false);
}
// scaling 
HEaaN::Ciphertext temp = HEaaN::Ciphertext(ctxt);

eval.mult(temp,0.01,temp);

approximate_derivated_ReLU(eval,btp,temp,ctxt_res);

if(ctxt_res.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
btp.bootstrap(ctxt_res,ctxt_res,false);
}

eval.mult(ctxt,ctxt_res,ctxt_res);
}

void approximate_Sigmoid(const HEaaN::HomEvaluator &eval,
const HEaaN::Bootstrapper &btp,
const HEaaN::Ciphertext &ctxt,
HEaaN::Ciphertext &ctxt_res,
const HEaaN::Real multiplier) {


std::vector<HEaaN::Real> coeffs = {
    5.000000000000016653e-01,
    2.494321294618147555e-01,
    1.197823256946907543e-16,
    -2.000661196427806765e-02,
    -1.270710679339945555e-16,
    1.697717910099649785e-03,
    2.630290438943759483e-17,
    -1.168492463579089727e-04,
    -2.811095069922330833e-18,
    6.146630386676855188e-06,
    1.845854028745479525e-19,
    -2.448328148258623576e-07,
    -8.141956968972177522e-21,
    7.464579518290908179e-09,
    2.561757863761087194e-22,
    -1.771673322298626663e-10,
    -6.003048484354804923e-24,
    3.330602094272319484e-12,
    1.081688587452919355e-25,
    -5.037886778776644381e-14,
    -1.535384417780968608e-27,
    6.215306306079636109e-16,
    1.748969521366997653e-29,
    -6.326251855024509795e-18,
    -1.622076451292634829e-31,
    5.363378911196642650e-20,
    1.238747833726888755e-33,
    -3.817019238997541022e-22,
    -7.858113203772922151e-36,
    2.294662555373150359e-24,
    4.168373269087639386e-38,
    -1.170928764520590443e-26,
    -1.857913505277351588e-40,
    5.089887674861345582e-29,
    6.980292112514033567e-43,
    -1.889162133147062535e-31,
    -2.214134629339008706e-45,
    5.993880945004973712e-34,
    5.929549721332120749e-48,
    -1.625493484092386132e-36,
    -1.338451808989476348e-50,
    3.762631819119272628e-39,
    2.537508942250889920e-53,
    -7.412840353683336507e-42,
    -4.017280369725787151e-56,
    1.237320048820821458e-44,
    5.265820953669235658e-59,
    -1.738200691480341513e-47,
    -5.645491827612264238e-62,
    2.036102057131350526e-50,
    4.864996745852294948e-65,
    -1.963304919403682157e-53,
    -3.285902886785858820e-68,
    1.530587735017802589e-56,
    1.674428696505288940e-71,
    -9.402381935211386360e-60,
    -6.050271092310076364e-75,
    4.379064812005916694e-63,
    1.381115621796050669e-78,
    -1.452691322600514629e-66,
    -1.497006196907696247e-82,
    3.057113553655033690e-70,
    0.000000000000000000e+00,
    -3.066553860586402329e-74,
};

const HEaaN::Math::ChebyshevCoefficients cheby_coeffs(coeffs, 16,
HEaaN::Math::PolynomialBasis::basic);

const HEaaN::Math::InputInterval input_interval(-25, 25);


HEaaN::Ciphertext tmp = HEaaN::Math::linearTransform(eval, btp, ctxt, input_interval);

    if (tmp.getLevel() < cheby_coeffs.level_cost + btp.getMinLevelForBootstrap()) {
        btp.bootstrap(tmp,tmp,false);
    }
    ctxt_res = HEaaN::Math::evaluateChebyshevExpansion(eval,btp,tmp,cheby_coeffs,1.0);

}