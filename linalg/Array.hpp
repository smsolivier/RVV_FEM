#pragma once 

#include "General.hpp"
#include <algorithm>

namespace fem 
{

/// Template class for wrapping std::vector. Mostly provides out of range checks 
template<class T=int> 
class Array {
public:
	/// default constructor 
	Array() { }
	/// initialize 
	/** \param N size of array 
	*/ 
	Array(int N) {_vector.resize(N); }
	/// initialize and set values 
	/** \param N size of array 
		\param val initial value 
	*/ 
	Array(int N, T val) {_vector.resize(N, val); } 
	/// construct from initializer list 
	Array(std::initializer_list<T> list) {
		for (auto i=list.begin(); i != list.end(); i++) {
			_vector.push_back(*i); 
		}
	}
	/// set from initializer list 
	void operator=(std::initializer_list<T> list) {
		Clear(); 
		for (auto i=list.begin(); i != list.end(); i++) {
			_vector.push_back(*i); 
		}
	}
	/// return the size of the array 
	int GetSize() const {return _vector.size(); }
	/// resize the array 
	void Resize(int N) {_vector.resize(N); }
	/// access to the array 
	T& operator[](int ind) {
		CHECKMSG(ind < _vector.size() && ind >= 0, 
			"ind = " << ind << ", size = " << _vector.size()); 
		return _vector[ind]; 
	}
	/// const access to the array 
	const T& operator[](int ind) const {
		CHECKMSG(ind < _vector.size() && ind >= 0, 
			"ind = " << ind << ", size = " << _vector.size()); 
		return _vector[ind]; 
	}
	/// set all elements to val 
	void operator=(T val) {
		for (int i=0; i<GetSize(); i++) {
			// (*this)[i] = val; 
			_vector[i] = val; 
		}
	}
	/// add to end of Array 
	void Append(T val) {_vector.push_back(val); }
	/// add an array to the back of this 
	void Append(const Array<T>& a) {
		int size = a.GetSize(); 
		for (int i=0; i<size; i++) {
			Append(a[i]); 
		}
	}
	/// clear contents of vector 
	void Clear() {_vector.clear(); }
	/// return the intersection of two arrays 
	void Intersection(const Array<T>& x, Array<T>& r) const {
		for (int i=0; i<x.GetSize(); i++) {
			const T& ival = x[i]; 
			for (int j=0; j<GetSize(); j++) {
				const T& jval = (*this)[j]; 
				if (jval == ival) {
					r.Append(ival); 
				}
			}
		}
	}
	/// remove repeated entries 
	void Unique() {
		std::sort(_vector.begin(), _vector.end()); 
		auto last = std::unique(_vector.begin(), _vector.end()); 
		_vector.erase(last, _vector.end()); 
	}
	/// sort from low to high 
	void Sort() {
		if (GetSize() > 0)
			std::sort(_vector.begin(), _vector.end()); 
	}
	/// check if an entry matches one in the array 
	bool In(T val) const {
		for (int i=0; i<GetSize(); i++) {
			if ((*this)[i] == val) return true; 
		}
		return false; 
	}
	/// test if two arrays are the same 
	bool operator==(const Array<T>& a) const {
		if (GetSize() != a.GetSize()) return false; 

		for (int i=0; i<GetSize(); i++) {
			if (a[i] != (*this)[i]) return false; 
		}
		return true; 
	}
	/// reverse order of array 
	void Transpose() {
		Array<T> tmp = (*this); 
		for (int i=0; i<tmp.GetSize(); i++) {
			(*this)[i] = tmp[tmp.GetSize() - 1 - i]; 
		}
	}
	/// print the Array
	void Print(std::ostream& out=std::cout) const {
		for (int i=0; i<GetSize(); i++) {
			out << (*this)[i] << " "; 
		}
		out << std::endl; 
	}
	/// access to the data 
	T* GetData() {return _vector.data(); }
	/// const access to the data 
	const T* GetData() const {return _vector.data(); }
private:
	/// vector that stores all the data 
	std::vector<T> _vector; 
}; 

template<> 
void Array<double>::Resize(int N); 
template<>
void Array<int>::Resize(int N); 

template<class T> 
std::ostream& operator<<(std::ostream& out, const Array<T>& a) {
	for (int i=0; i<a.GetSize(); i++) {
		out << a[i] << " "; 
	}
	return out; 
}

} // end namespace fem 