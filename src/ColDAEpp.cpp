// ColDAEpp.cpp: Definiert den Einstiegspunkt für die Anwendung.
//

#include "ColDAEpp.hpp"

using namespace std;





#include <vector>


#include <iostream>

void test(int& sdf) {

}
int main()
{
	arrData1<int> ad = { 3,4,56,6 };

	arrRef1<int> ar = ad;
	for (auto a : ar)
		std::cout << a << " ";
	std::cout << std::endl;
	arrRef1<int> ar2 = ad.sub(3);
	for (auto a : ar2)
		std::cout << a << " " ;
	std::cout << std::endl;

	test(ar(1));

	ar(1) = 567;
	ad(1) = 567;

	for (auto a : ar)
		std::cout << a << " ";

	std::cin.get();
	return 0;
}
