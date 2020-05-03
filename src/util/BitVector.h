#pragma once

#include <vector>

namespace masc {
	namespace unfolding {

		using namespace std;

		//a long stream of bits
		class BitVector {
		public:

			BitVector(const vector<bool>& values);
			BitVector(int size);//initialize all to 0s

			int size() const { return _N; }
			void on(int i); //turn on i-th bit
			void off(int i); //turn off i-th bit

			//operator overloading
			bool operator[](int i) const;
			bool operator<(const BitVector& other) const; //comparator
			friend ostream& operator<<(ostream& out, const BitVector& v);
			bool operator==(const BitVector& other) const;

			//convert BitVector to a vector of bool
			void tovector(vector<bool>& values) const;

		protected:
			vector<unsigned int> data;
			const static int _S; //this is sizeof(unsigned int)
			int _N; //number of bits stored
			friend struct hash_BitVector;
		};

		struct hash_BitVector {
			//for hashing
			size_t operator()(const BitVector &v) const
			{
				size_t code = 0x2D2816FE;
				for (int i : v.data) return code * 31 + i; //hash<int>()(i);
				return code;
			}
		};
    }
}
