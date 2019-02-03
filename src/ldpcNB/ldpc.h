#pragma once

#include <iostream>
#include <exception>
#include <string>

using namespace std;

namespace ldpc {
    class GFMat;
    class LdpcCode;

    void rankGFMatrix(GFMat* mat);
}

using namespace ldpc;
