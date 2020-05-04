#include "BitVector.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace masc {
	namespace unfolding {

		const int BitVector::_S = sizeof(int) * 8; //number of bits for int

		BitVector::BitVector(const vector<bool>& values)
			:_N((int)values.size())
		{
			data = vector<unsigned int>((int)ceil(_N*1.0f / _S), 0);

			for (int i = 0; i < _N; i++) {
				//cout<<"trying i="<<i<<endl;
				if (values[i]) on(i);
				// if((*this)[129]){
				//   cout<<"129 is on why? "<<i<<endl;
				//   exit(1);
				// }
			}
		}

		BitVector::BitVector(int size)
			:_N(size)
		{
			data = vector<unsigned int>((int)ceil(_N*1.0f / _S), 0);
		}

		void BitVector::on(int i) {
			assert(i < _N);
			data[i / _S] |= (1 << (i%_S));
		}

		void BitVector::off(int i) {
			assert(i < _N);
			data[i / _S] &= ~(1 << (i%_S));
		}

		bool BitVector::operator[](int i) const {
			assert(i < _N);
			return data[i / _S] & (1 << (i%_S));
		}

		bool BitVector::operator<(const BitVector& other) const {
			int size = data.size();
			for (int i = size - 1; i >= 0; i--) {
				if (data[i] != other.data[i]) return data[i] < other.data[i];
			}
			return false; //they are the same
		}

		ostream& operator<<(ostream& out, const BitVector& v) {
			out << "(";
			for (int i = 0; i < v._N; i++) out << v[i] << ((i < v._N - 1) ? "," : ")");
			return out;
		}

		bool BitVector::operator==(const BitVector& other) const {
			if (other._N != _N) return false;
			int size = data.size();
			for (int i = 0; i < size; i++) if (data[i] != other.data[i]) return false;
			return true;
		}

		void BitVector::tovector(vector<bool>& values) const {
			values.resize(_N); //ensure that there is enough space
			for (int i = 0; i < _N; i++) values[i] = (*this)[i];
		}

    }
}
