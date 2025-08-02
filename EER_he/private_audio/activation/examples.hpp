////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2024 CryptoLab, Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of CryptoLab, Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>

#include "HEaaN/HEaaN.hpp"


inline double randNum() {
    static std::default_random_engine gen{std::random_device()()};
    std::uniform_real_distribution<double> dist(-1.0L, 1.0L);
    return dist(gen);
}

inline void fillRandomComplex(HEaaN::Message &msg) {
    for (size_t i = 0; i < msg.getSize(); ++i) {
        msg[i].real(randNum());
        msg[i].imag(randNum());
    }
}

inline void fillRandomReal(HEaaN::Message &msg, std::optional<size_t> num) {
    size_t length = num.has_value() ? num.value() : msg.getSize();
    size_t idx = 0;
    for (; idx < length; ++idx) {
        //msg[idx].real(0.5);
        msg[idx].real(randNum()*25);
        msg[idx].imag(0.0);
    }
    // If num is less than the size of msg,
    // all remaining slots are zero.
    for (; idx < msg.getSize(); ++idx) {
        msg[idx].real(0.0);
        msg[idx].imag(0.0);
    }
}

inline void printMessage(const HEaaN::Message &msg, bool is_complex = true,
                         size_t start_num = 2, size_t end_num = 2,
                         const int precision = 10) {
    const size_t msg_size = msg.getSize();
    std::cout.precision(precision);
    std::cout << "[ ";
    for (size_t i = 0; i < start_num; ++i) {
        if (is_complex)
            std::cout << msg[i] << ", ";
        else
            std::cout << msg[i].real() << ", ";
    }
    std::cout << "..., ";
    for (size_t i = end_num; i > 1; --i) {
        if (is_complex)
            std::cout << msg[msg_size - i] << ", ";
        else
            std::cout << msg[msg_size - i].real() << ", ";
    }
    if (is_complex)
        std::cout << msg[msg_size - 1] << " ]" << std::endl;
    else
        std::cout << msg[msg_size - 1].real() << " ]" << std::endl;
}

inline std::string presetNamer(const HEaaN::ParameterPreset preset) {
    switch (preset) {
    case HEaaN::ParameterPreset::FVa:
        return "FVa";
    case HEaaN::ParameterPreset::FVb:
        return "FVb";
    case HEaaN::ParameterPreset::FGa:
        return "FGa";
    case HEaaN::ParameterPreset::FGb:
        return "FGb";
    case HEaaN::ParameterPreset::FTa:
        return "FTa";
    case HEaaN::ParameterPreset::FTb:
        return "FTb";
    case HEaaN::ParameterPreset::ST19:
        return "ST19";
    case HEaaN::ParameterPreset::ST14:
        return "ST14";
    case HEaaN::ParameterPreset::ST11:
        return "ST11";
    case HEaaN::ParameterPreset::ST8:
        return "ST8";
    case HEaaN::ParameterPreset::ST7:
        return "ST7";
    case HEaaN::ParameterPreset::SS7:
        return "SS7";
    case HEaaN::ParameterPreset::SD3:
        return "SD3";
    case HEaaN::ParameterPreset::CUSTOM:
        return "CUSTOM";
    case HEaaN::ParameterPreset::FVc:
        return "FVc";
    case HEaaN::ParameterPreset::FGd:
        return "FGd";
    case HEaaN::ParameterPreset::SGd0:
        return "SGd0";
    default:
        throw std::invalid_argument("Not supported parameter");
    }
}

void printCopyright() {
    if (std::filesystem::exists("copyright")) {
        std::ifstream copyright_info;
        copyright_info.open("copyright");

        while (!copyright_info.eof()) {
            std::string line;
            std::getline(copyright_info, line);
            std::cout << line << std::endl;
        }
    } else {
        std::cout << "ⓒ 2024 CryptoLab, Inc. All rights reserved.\n"
                  << std::endl;
    }
}
