#include <iostream>
#include "lib/string.h"
#include "lib/list.h"
#include "lib/unordered_map.h"
#include "lib/deque.h"
#include "lib/geometry.h"
#include "lib/biginteger.h"


int main() {
    BigInteger kek = (long long)1e18;
    BigInteger kek2 = (long long)1e18;
    BigInteger kek3 = kek * kek2;
    kek3 *= kek3;
    kek3 *= kek3;
    kek3 *= kek3;
    kek3 *= kek3;
    std::cout << kek3 << '\n';
    return 0;
}
