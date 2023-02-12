#include <iostream>

template<int n, typename I>
auto align_up(const I i) -> I {
	return i%n ? i/n*n+n : i;
}

int main(){
	for(int i=0; i<9; i++)
		std::cout << align_up<4>(i) << std::endl;
}
