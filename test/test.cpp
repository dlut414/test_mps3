// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>

template <unsigned n>
class factorial {
public:
	enum { value = n* factorial<n - 1>::value, };
};

template <>
class factorial < 0 > {
public:
	enum { value = 1, };
};

template <unsigned n, unsigned m>
class H {
public:
	enum { value = factorial<n + m - 1>::value / (factorial<m>::value * factorial<n - 1>::value), };
};

template <typename R, unsigned D, unsigned P>
class Polynomial_A_ {
public:
	enum { value = H<D, P>::value + Polynomial_A_<R, D, P - 1>::value, };
	template <unsigned N>	static __forceinline R Power(const R& x) { return x*Power<N - 1>(x); }
	template <>				static __forceinline R Power<0>(const R& x) { return 1; }
};
template <typename R, unsigned D>	class Polynomial_A_<R, D, 0> { public: enum { value = 0, }; };

template <typename R, unsigned D, unsigned P>	class Polynomial_A : public Polynomial_A_<R, D, P> {};

template <typename R, unsigned P>				class Polynomial_A<R, 1, P> : public Polynomial_A_<R, 1, P>{
private:
	template <unsigned P_ = 1>	static __forceinline void Gen_(const R& x, R* const out) { out[P_ - 1] = Power<P_>(x); Gen_<P_ + 1>(x, out); }
	template <>					static __forceinline void Gen_<P>(const R& x, R* const out) { out[P - 1] = Power<P>(x); }
	template <unsigned P_ = 1>	static __forceinline void Gen_(const R& s, const R& x, R* const out) { out[P_ - 1] = Power<P_>(x / s); Gen_<P_ + 1>(s, x, out); }
	template <>					static __forceinline void Gen_<P>(const R& s, const R& x, R* const out) { out[P - 1] = Power<P>(x / s); }
public:
	static __forceinline void Gen(const R& x, R* const out) { Gen_(x, out); }
	static __forceinline void Gen(const R& s, const R& x, R* const out) { Gen_(s, x, out); }
};

template <typename R, unsigned P>				class Polynomial_A<R, 2, P> : public Polynomial_A_<R, 2, P>{
private:
	template <unsigned P_ = 1, unsigned PX = P_, unsigned I = 0>
	struct Gen__ {
		static __forceinline void Gen(const R* const x, R* const out) { out[I] = Power<PX>(x[0])*Power<P_ - PX>(x[1]); Gen__<P_, PX - 1, I + 1>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out) { out[I] = Power<PX>(x[0] / s)*Power<P_ - PX>(x[1] / s); Gen__<P_, PX - 1, I + 1>::Gen(s, x, out); }
	};
	template <unsigned P_, unsigned I>
	struct Gen__<P_, 0, I> {
		static __forceinline void Gen(const R* const x, R* const out) { out[I] = Power<P_>(x[1]); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out) { out[I] = Power<P_>(x[1] / s); }
	};
	template <unsigned P_ = 1, unsigned I = 0, bool Over = (P==P_)>
	struct Gen_ {
		static __forceinline void Gen(const R* const x, R* const out) { Gen__<P_, P_, I>::Gen(x, out); Gen_<P_ + 1, I + H<2, P_>::value>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out) { Gen__<P_, P_, I>::Gen(s, x, out); Gen_<P_ + 1, I + H<2, P_>::value>::Gen(s, x, out); }
	};
	template <unsigned P_, unsigned I>
	struct Gen_<P_, I, true> {
		static __forceinline void Gen(const R* const x, R* const out) { Gen__<P_, P_, I>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out) { Gen__<P_, P_, I>::Gen(s, x, out); }
	};
public:
	static __forceinline void Gen(const R* const x, R* const out) { Gen_<>::Gen(x, out); }
	static __forceinline void Gen(const R& s, const R* const x, R* const out) { Gen_<>::Gen(s, x, out); }
};

template <typename R, unsigned D, unsigned P>
class Derivative_A_ {
public:
	template <unsigned N>	static __forceinline R Power(const R& x) { return x*Power<N - 1>(x); }
	template <>				static __forceinline R Power<0>(const R& x) { return 1; }
};

template <typename R, unsigned D, unsigned P>	class Derivative_A : public Derivative_A_<R, D, P> {};

template <typename R, unsigned P>				class Derivative_A<R, 1, P> : public Derivative_A_<R, 1, P>{
private:
	template <unsigned Der, unsigned P_ = 1, bool Zero = (Der>P_) >
	struct Gen_ {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = Der*Power<P_ - Der>(x); Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = Der*Power<P_ - Der>(x / s)*Power<Der>(1 / s); Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(s, x, out); }
	};
	template <unsigned Der, unsigned P_>
	struct Gen_<Der, P_, true> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = 0; Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = 0; Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(s, x, out); }
	};
	template <unsigned Der>
	struct Gen_<Der, P, true> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = 0; }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = 0; }
	};
	template <unsigned Der>
	struct Gen_<Der, P, false> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = Der*Power<P_ - Der>(x); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = Der*Power<P_ - Der>(x / s)*Power<Der>(1 / s); }
	};
public:
	template <unsigned Der> static __forceinline void Gen(const R& x, R* const out) { Gen_<Der>::Gen(x, out); }
	template <unsigned Der> static __forceinline void Gen(const R& s, const R& x, R* const out) { Gen_<Der>::Gen(s, x, out); }
};

template <unsigned N>
class TEST {
public:
	template <unsigned P, unsigned I, bool Specialize = (N==P)>		struct test					{ static __forceinline void Run() { std::cout << 0 << std::endl; } };
	template <unsigned P, unsigned I>								struct test <P, I, true>	{ static __forceinline void Run() { std::cout << 1 << std::endl; } };
	template <unsigned P>											struct test <P, 1, true>	{ static __forceinline void Run() { std::cout << 0 << std::endl; } };
};

int _tmain(int argc, _TCHAR* argv[])
{
	TEST<12>::test<12, 2>::Run();

	//typedef Polynomial_A<float, 2, 4> Poly;

	//typedef Eigen::Matrix<float, 2, 1> vec21;
	//typedef Eigen::Matrix<float, Poly::value, 1> vecPoly;

	////vec21 in;
	//vecPoly out;
	//in << 1, 2;

	////Polynomial_A<float, 1, 2>::Gen(3, out.data());
	////std::cout << out << std::endl;

	//Poly::Gen(in.data(), out.data());
	//std::cout << " input is: " << std::endl; 
	//std::cout << in << std::endl;
	//std::cout << std::endl;
	//std::cout << " output contains: " << Poly::value << " values " << std::endl;
	//std::cout << out << std::endl;

	std::cout << " Derivative_A test: " << std::endl;
	typedef Polynomial_A<float, 1, 1> Poly;
	typedef Eigen::Matrix<float, Poly::value, 1> vecPoly;
	typedef Derivative_A<float, 1, 1> Der;
	float derIn = 2;
	vecPoly out;
	Der::Gen<1>(derIn, out.data());
	std::cout << out << std::endl;

	return 0;
}

