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
	template <unsigned N, unsigned M>	static __forceinline R A() { return (factorial<N>::value)/(factorial<N-M>::value); }
};

template <typename R, unsigned D, unsigned P>	class Derivative_A : public Derivative_A_<R, D, P> {};

template <typename R, unsigned P>				class Derivative_A<R, 1, P> : public Derivative_A_<R, 1, P>{
private:
	template <unsigned Der, unsigned P_ = 1, bool Zero = (Der>P_), bool Over = (P_==P)>
	struct Gen_ {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = A<P_,Der>()*Power<P_ - Der>(x); Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = A<P_, Der>()*Power<P_ - Der>(x / s)*Power<Der>(1 / s); Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(s, x, out); }
	};
	template <unsigned Der, unsigned P_, bool Over>
	struct Gen_<Der, P_, true, Over> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = 0; Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = 0; Gen_<Der, P_ + 1, (Der>(P_ + 1))>::Gen(s, x, out); }
	};
	template <unsigned Der, unsigned P_>
	struct Gen_<Der, P_, true, true> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = 0; }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = 0; }
	};
	template <unsigned Der, unsigned P_>
	struct Gen_<Der, P_, false, true> {
		static __forceinline void Gen(const R& x, R* const out) { out[P_ - 1] = A<P_, Der>()*Power<P_ - Der>(x); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_ - 1] = A<P_, Der>()*Power<P_ - Der>(x / s)*Power<Der>(1 / s); }
	};
public:
	template <unsigned Der> static __forceinline void Gen(const R& x, R* const out) { Gen_<Der>::Gen(x, out); }
	template <unsigned Der> static __forceinline void Gen(const R& s, const R& x, R* const out) { Gen_<Der>::Gen(s, x, out); }
};

	template <typename R, unsigned P>				class Derivative_A<R,2,P> : public Derivative_A_<R,2,P> {
	private:
		template <unsigned DerX, unsigned DerY, unsigned P_ = 1, unsigned PX = P_, unsigned I = 0, bool Zero = (DerX>PX||DerY>(P_-PX))>
		struct Gen__									{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = A<PX,DerX>()*A<P_-PX,DerY>()*Power<PX-DerX>(x[0])*Power<P_-PX-DerY>(x[1]); Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = A<PX,DerX>()*A<P_-PX,DerY>()*Power<DerX+DerY>(1/s)*Power<PX-DerX>(x[0]/s)*Power<P_-PX-DerY>(x[1]/s); Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned PX, unsigned I>
		struct Gen__ <DerX,DerY,P_,PX,I,true>			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = 0; Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = 0; Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen__<DerX,DerY,P_,0,I,true> 			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = 0; }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = 0; }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen__<DerX,DerY,P_,0,I,false> 			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = A<0,DerX>()*A<P_-0,DerY>()*Power<0-DerX>(x[0])*Power<P_-0-DerY>(x[1]); }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = A<0,DerX>()*A<P_-0,DerY>()*Power<DerX+DerY>(1/s)*1*Power<P_-0-DerY>(x[1]/s); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_ = 1, unsigned I = 0, bool Over = (P_==P)>
		struct Gen_							{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<DerX,DerY,P_,P_,I>::Gen(x, out); Gen_<DerX,DerY,P_+1,I+H<2,P_>::value>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<DerX,DerY,P_,P_,I>::Gen(s, x, out); Gen_<DerX,DerY,P_+1,I+H<2,P_>::value>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen_<DerX,DerY,P_,I,true>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<DerX,DerY,P_,P_,I>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<DerX,DerY,P_,P_,I>::Gen(s, x, out); }
		};
	public:
		template <unsigned DerX, unsigned DerY>	static __forceinline void Gen(const R* const x, R* const out) { Gen_<DerX,DerY>::Gen(x, out); }
		template <unsigned DerX, unsigned DerY>	static __forceinline void Gen(const R& s, const R* const x, R* const out) { Gen_<DerX,DerY>::Gen(s, x, out); }
	};

template <unsigned N>
class TEST {
public:
	template <unsigned P, unsigned I, bool Specialize = (N==P)>		struct test					{ static __forceinline void Run() { std::cout << 0 << std::endl; } };
	template <unsigned P, unsigned I>								struct test <P, I, true>	{ static __forceinline void Run() { std::cout << 1 << std::endl; } };
	template <unsigned P>											struct test <P, 1, true>	{ static __forceinline void Run() { std::cout << 0 << std::endl; } };
};

template <typename R, unsigned D>
class LinkCell {
	typedef Eigen::Matrix<R, D, 1>		vecD;
	typedef Eigen::Matrix<int, D, 1>	iVecD;
public:
	template <unsigned D_ = 0, bool Over = (D_ == (D - 1))> struct Convert {
		template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out) { out[D_] = static_cast<V>(in[D_]); Convert<D_ + 1>::Gen(in, out); }
	};
	template <unsigned D_>								struct Convert < D_, true > {
		template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out) { out[D_] = static_cast<V>(in[D_]); }
	};

	static __forceinline const iVecD iCoord(const vecD& p) { iVecD ret; Convert<>::Gen(p.data(), ret.data()); return ret; }
};

int _tmain(int argc, _TCHAR* argv[])
{
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

	std::cout << " Derivative_A 1D test: " << std::endl;
	typedef Polynomial_A<float, 1, 4> Poly1D;
	typedef Eigen::Matrix<float, Poly1D::value, 1> vecPoly1D;
	typedef Derivative_A<float, 1, 4> Der1D;
	float derIn1D = 2;
	vecPoly1D out1D;
	Der1D::Gen<3>(2, derIn1D, out1D.data());
	std::cout << out1D << std::endl;

	std::cout << " Derivative_A 2D test: " << std::endl;
	typedef Polynomial_A<float, 2, 4> Poly2D;
	typedef Eigen::Matrix<float, Poly2D::value, 1> vecPoly2D;
	typedef Derivative_A<float, 2, 4> Der2D;
	Eigen::Matrix<float, 2, 1> derIn2D;
	derIn2D << 1, 2;
	vecPoly2D out;
	Der2D::Gen<1,0>(2, derIn2D.data(), out.data());
	std::cout << out << std::endl;

	std::cout << " converter test: " << std::endl;
	Eigen::Matrix<float, 2, 1> input;
	input << 1.1f, 2.1f;
	std::cout << LinkCell<float, 2>::iCoord(input) << std::endl;

	std::cout << " test over " << std::endl;

	return 0;
}

