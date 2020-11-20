#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

int main() {
    ZZ a = ZZ(14);
    ZZ b = ZZ(15); 

    cout << (a + b) << "\n";
}