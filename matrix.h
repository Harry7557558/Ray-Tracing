#pragma once
#pragma warning(disable: 4244)	// conversion from 'type1' to 'type2', possible loss of data
#pragma warning(disable: 4018)	// signed/unsigned mismatch, note that sometimes cause problems


#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <queue>
#include <stack>

#include <cstdlib>
#include <ctime>

using namespace std;

#ifndef _INC_FRACTION

#define _INC_FRACTION


static unsigned _Myt_Fraction_GCF(unsigned a, unsigned b)
{
	if (a == 0 || b == 0) return 0;
	unsigned c = -1;
	while (c != 0) c = a % b, a = b, b = c;
	return a;
}
static unsigned _Myt_Fraction_LCM(unsigned a, unsigned b)
{
	if (a == 0 || b == 0) return 0;
	return b / _Myt_Fraction_GCF(a, b) * a;
}
static unsigned _Myt_Fraction_String_to_Num(string a) {
	unsigned n = 0;
	while (a.size() != 0) {
		n *= 10, n += a[0] - 48;
		a.erase(0, 1);
	}
	return n;
}


class fraction {
	unsigned den;
	int num;
	inline friend void _Myt_Fraction_RFCD(fraction &a, fraction &b)
	{
		unsigned c = _Myt_Fraction_LCM(a.den, b.den);
		a.num = c / a.den*a.num;
		b.num = c / b.den*b.num;
		a.den = b.den = c;
	}
	static bool _Myt_Fraction_IssDec(fraction a) {
		int q = log10(abs(a.num)) + log10(a.den) + 2;
		unsigned char s = 0, t = 0;
		while (!(a.den & 1)) a.den >>= 1, t++;
		while (!(a.den % 5)) a.den /= 5, s++;
		if (a.den != 1) return 0;
		t = t > s ? t : s;
		if (t <= 2 && log10(abs(a.num)) < 3) return 1;
		q -= t;
		q -= log10(abs(a.num) / a.den);
		return q >= 0;
	}
public:
	fraction() :num(0), den(1) {}
	fraction(int numerator, int denominator) {
		if (denominator < 0) numerator = -numerator, denominator = -denominator;
		num = numerator, den = denominator;
		simplify();
	}
	fraction(const int &a) {
		num = a, den = 1;
	}
	fraction(const double &a) {

	}
	fraction(const string &a) {
		*this = a;
	}
	fraction& operator = (const int &a) {
		num = a, den = 1;
		return *this;
	}
	fraction& operator = (const double &a) {
		fraction f;
		int w = a;
		if (a == w) {
			f.num = w, f.den = 1;
			return f;
		}
		unsigned n = 1, d = 2;
		double s;
		do {
			s = (w*d + n) / d;
			if (s > a) {
				d *= 2, n = n * 2 - 1;

			}
			if (s < a) {
				d *= 2, n = n * 2 + 1;
			}
		} while (s != a);
	}
	void operator = (string a) {
		bool sign = 0;
		if (a[0] == '-') sign = 1, a.erase(0, 1);
		if (a[a.find('/', 0) + 1] == '-' && a.find('/', 0) != 0) sign ^= 1;
		if (a.find('/', 0) == -1) {
			if (a.find('.', 0) != -1) {
				for (unsigned i = 0; i < a.size(); i++) {
					if ((a[i] < 47 || a[i] > 57) && a[i] != '.') a.erase(i, 1);
				}
				for (unsigned i = a.find('.', 0) + 1; i < a.size(); i++) {
					if (a[i] == '.') a.erase(i, 1), i--;
				}
				den = 1;
				for (unsigned i = a.size() - a.find('.', 0) - 1; i > 0; i--) {
					den *= 10;
				}
				a.erase(a.find('.', 0), 1);
				num = _Myt_Fraction_String_to_Num(a);
				if (sign) num = -num;
				simplify();
				return;
			}
			for (unsigned i = 0; i < a.size(); i++) {
				if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
			}
			num = _Myt_Fraction_String_to_Num(a), den = 1;
			if (sign) num = -num;
			simplify();
			return;
		}
		for (unsigned i = 0; i < a.size(); i++) {
			if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
		}
		num = _Myt_Fraction_String_to_Num(a.substr(0, a.find('/', 0)));
		den = _Myt_Fraction_String_to_Num(a.substr(a.find('/', 0) + 1, a.size() - a.find('/', 0)));
		if (sign) num = -num;
		simplify();
	}
	inline friend ostream& operator << (ostream& os, fraction a)
	{
		if (a.den == 0) {
			if (a.num > 0) os << "#INF";
			else if (a.num < 0) os << "#-INF";
			else if (a.num == 0) os << "#NAF";
			return os;
			// "#INF" will become "#NAN" after any calculation. 
		}
		if (a.num == 0) {
			os << "0";
			return os;
		}
		if (a.den == 1) {
			os << a.num;
			return os;
		}
		a.simplify();
		/*if (_Myt_Fraction_IssDec(a)) {
			os << (double)a.num / (double)a.den;
			return os;
		}*/
		os << a.num << "/" << a.den;
		return os;
	}
	inline friend istream& operator >> (istream& is, fraction &a)
	{
		string n; is >> n;
		a = n;
		return is;
	}
	inline void simplify()
	{
		if (den == 0) {
			if (num > 0) num = 1;
			if (num < 0) num = -1;
			return;
		}
		if (num == 0) {
			den = 1; return;
		}
		int m = _Myt_Fraction_GCF(abs(num), den);
		if (m == 1) return;
		den /= m; num /= m;
	}
	inline bool integer() {
		return den == 1;
	}
	inline unsigned denominator() {
		return den;
	}
	inline int numerator() {
		return num;
	}
	inline friend bool operator == (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		a.simplify(), b.simplify();
		return (a.den == b.den && a.num == b.num);
	}
	inline friend bool operator != (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		return !(a == b);
	}
	inline friend bool operator == (fraction a, int b)
	{
		return (b*a.den == a.num);
	}
	inline friend bool operator != (fraction a, int b)
	{
		return (b*a.den != a.num);
	}
	inline friend bool operator < (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		_Myt_Fraction_RFCD(a, b);
		return (a.num < b.num);
	}
	inline friend bool operator > (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		_Myt_Fraction_RFCD(a, b);
		return (a.num > b.num);
	}
	inline friend bool operator <= (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		_Myt_Fraction_RFCD(a, b);
		return (a.num <= b.num);
	}
	inline friend bool operator >= (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		_Myt_Fraction_RFCD(a, b);
		return (a.num >= b.num);
	}
	inline friend bool operator > (fraction a, int b) {
		return b * int(a.den) < a.num;
	}
	inline friend bool operator < (fraction a, int b) {
		return b * int(a.den) > a.num;
	}
	inline friend bool operator >= (fraction a, int b) {
		return b * a.den <= a.num;
	}
	inline friend bool operator <= (fraction a, int b) {
		return b * a.den >= a.num;
	}
	inline friend bool operator == (fraction a, double b) {
		return double(a.num) / a.den == b;
	}
	inline friend bool operator != (fraction a, double b) {
		return double(a.num) / a.den != b;
	}
	inline friend bool operator > (fraction a, double b) {
		return double(a.num) / a.den > b;
	}
	inline friend bool operator >= (fraction a, double b) {
		return double(a.num) / a.den >= b;
	}
	inline friend bool operator < (fraction a, double b) {
		return double(a.num) / a.den < b;
	}
	inline friend bool operator <= (fraction a, double b) {
		return double(a.num) / a.den <= b;
	}
	inline friend bool operator ! (fraction a) {
		return !a.num;
	}
	inline friend fraction operator + (fraction a, int b)
	{
		a.num += b * a.den;
		a.simplify();
		return a;
	}
	inline friend void operator += (fraction &a, int b)
	{
		a.num += b * a.den;
		a.simplify();
	}
	inline friend fraction operator + (int a, fraction b)
	{
		b.num += a * b.den;
		b.simplify();
		return b;
	}
	inline friend fraction operator + (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) {
			a.den = a.num = 0; return a;
		}
		if (a.num == 0) {
			b.simplify(); return b;
		}
		if (b.num == 0) {
			a.simplify(); return a;
		}
		_Myt_Fraction_RFCD(a, b);
		a.num += b.num;
		a.simplify();
		return a;
	}
	inline friend void operator += (fraction &a, fraction b)
	{
		if (a.den == 0 || b.den == 0) {
			a.den = a.num = 0; return;
		}
		if (b.num == 0) {
			a.simplify(); return;
		}
		if (a.num == 0) {
			b.simplify(); a = b; return;
		}
		_Myt_Fraction_RFCD(a, b);
		a.num += b.num;
		a.simplify();
	}
	inline friend fraction& operator ++ (fraction &f) {
		f.num += f.den;
		return f;
	}
	inline friend fraction operator ++ (fraction &f, int) {
		f.num += f.den;
		return f;
	}
	inline friend fraction operator - (fraction a)
	{
		a.num = -a.num;
		a.simplify();
		return a;
	}
	inline friend fraction operator - (fraction a, int b)
	{
		a.num -= b * a.den;
		a.simplify();
		return a;
	}
	inline friend fraction operator - (int a, fraction b)
	{
		b.num -= a * b.den;
		b.simplify();
		return b;
	}
	inline friend void operator -= (fraction &a, int b)
	{
		a.num -= b * a.den;
		a.simplify();
	}
	inline friend fraction operator - (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) {
			a.den = a.num = 0; return a;
		}
		if (a.num == 0) {
			return -b;	// Already simplified in "operator -" 
		}
		if (b.num == 0) {
			a.simplify(); return a;
		}
		_Myt_Fraction_RFCD(a, b);
		a.num -= b.num;
		a.simplify();
		return a;
	}
	inline friend void operator -= (fraction &a, fraction b)
	{
		if (a.den == 0 || b.den == 0) {
			a.den = a.num = 0; return;
		}
		if (b.num == 0) {
			a.simplify(); return;
		}
		if (a.num == 0) {
			a = -b; return;
		}
		_Myt_Fraction_RFCD(a, b);
		a.num -= b.num;
		a.simplify();
	}
	inline friend fraction& operator -- (fraction &f) {
		f.num -= f.den;
		return f;
	}
	inline friend fraction operator -- (fraction &f, int) {
		f.num -= f.den;
		return f;
	}
	inline friend fraction operator * (fraction a, int b)
	{
		a.num *= b;
		a.simplify();
		return a;
	}
	inline friend fraction operator * (int a, fraction b)
	{
		b.num *= a;
		b.simplify();
		return b;
	}
	inline friend void operator *= (fraction &a, int b)
	{
		a.num *= b;
		a.simplify();
	}
	inline friend fraction operator * (fraction a, fraction b)
	{
		a.den *= b.den, a.num *= b.num;
		a.simplify();
		return a;
	}
	inline friend void operator *= (fraction &a, fraction b)
	{
		a.den *= b.den, a.num *= b.num;
		a.simplify();
	}
	inline friend fraction operator / (fraction a, int b)
	{
		if (b < 0) a.num = -a.num, b = -b;
		a.den *= b;
		a.simplify();
		return a;
	}
	inline friend void operator /= (fraction &a, int b)
	{
		if (b < 0) a.num = -a.num, b = -b;
		a.den *= b;
		a.simplify();
	}
	inline friend fraction operator / (int a, fraction b)
	{
		bool sign = a < 0;
		if (sign) a = -a;
		if (b.num < 0) b.num = -b.num, sign ^= 1;
		b.den *= a;
		int k = b.den;
		if (sign) k = -k;
		b.den = b.num, b.num = k;
		b.simplify();
		return b;
	}
	inline friend fraction operator / (fraction a, fraction b)
	{
		a.num *= b.den, a.den *= abs(b.num);
		if (b.num < 0) a.num = -a.num;
		a.simplify();
		return a;
	}
	inline friend void operator /= (fraction &a, fraction b)
	{
		a.num *= b.den, a.den *= abs(b.num);
		if (b.num < 0) a.num = -a.num;
		a.simplify();
	}
	inline friend fraction operator ^ (fraction a, int k)
	{
		if (k == 0) {
			a.num = a.den = 1;
			return a;
		}
		a.simplify();
		if (a.num == 0) return a;
		if (k < 0) {
			a = ~a;
			k = -k;
		}
		fraction b; b.den = b.num = 1;
		for (int i = 0; i < k; i++) {
			b.den *= a.den;
			b.num *= a.num;
		}
		return b;
	}
	inline friend fraction pow(fraction a, int k) {
		return a ^ k;
	}
	inline friend void operator ^= (fraction &a, unsigned k)
	{
		a = a ^ k;
	}
	inline friend fraction operator ~ (fraction a)
	{
		bool s = a.num < 0;
		fraction b;
		b.num = a.den, b.den = abs(a.num);
		if (s) b.num = -b.num;
		b.simplify();
		return b;
	}
};


#endif

#define TeX_endl "\\\\[8pt]\n";
#define TeX_rightarrow "\\Longrightarrow"

#ifndef _INC_MATRIX

#define _INC_MATRIX

namespace Matrix_type {
	struct _M_type {
		unsigned char type = 0;
	};
	inline _M_type operator & (_M_type a, _M_type b) {
		a.type |= b.type;
		return a;
	}
	inline _M_type operator | (_M_type a, _M_type b) {
		a.type |= b.type;
		return a;
	}
	inline _M_type operator + (_M_type a, _M_type b) {
		a.type |= b.type;
		return a;
	}
	inline _M_type M_type(unsigned char a) {
		_M_type n; n.type = a; return n;
	}

#define _M_Identity M_type(0x01)
#define _M_Row M_type(0x02)
#define _M_Column M_type(0x04)

}

#include <initializer_list>
typedef unsigned char byte;

template<typename T> class matrix;

template<typename T> class Vector {
protected:
	unsigned l;
	T* data;
	friend class matrix<T>;
public:
	inline Vector() : l(0), data(0) {}
	Vector(const int &length) {
		l = length;
		data = new T[length];
	}
	Vector(const initializer_list<T> &q) {
		l = q.size();
		data = new T[l];
		for (int i = 0; i < l; i++) data[i] = (q.begin())[i];
	}
	Vector(const Vector<T> &V) {
		l = V.l; data = new T[l];
		for (unsigned i = 0; i < l; i++) data[i] = (V.data)[i];
	}
	~Vector() {
		l = 0;
		if (data != 0) delete data;
	}
	Vector<T>& operator = (const Vector<T> &V) {
		if (l != V.l) {
			l = V.l;
			if (data != 0) delete data, data = new T[l];
		}
		if (data == 0) data = new T[l];
		for (int i = 0; i < l; i++) data[i] = (V.data)[i];
		return *this;
	}
	Vector<T>& operator = (const initializer_list<T> &V) {
		if (l != V.size()) {
			l = V.size();
			if (data != 0) delete data, data = new T[l];
		}
		if (data == 0) data = new T[l];
		for (int i = 0; i < l; i++) data[i] = (V.begin())[i];
		return *this;
	}
	inline T& operator [] (const unsigned &k) {
		return data[k];
	}
	inline unsigned size() const { return l; }
	inline unsigned length() const { return l; }
	inline unsigned dimension() const { return l; }

	friend ostream& operator << (ostream& os, const Vector<T> &V)
	{
		os << "\\begin{bmatrix}";
		for (unsigned i = 0; i < V.length(); i++) {
			os << V.data[i] << "\\\\[3pt]";
		}
		os << "\\end{bmatrix}";
		return os;
	}

	//friend void solve_q(matrix<T> A, Vector<T> &X, const Vector<T> &Y);
};

template<typename T> class matrix {
protected:
	unsigned h, w;
	T* data;
public:
	inline matrix() :h(0), w(0), data(0) {}
	matrix(int length) {
		h = w = length;
		if (length == 0) {
			data = 0; return;
		}
		data = new T[length*length];
		for (int i = length * length; i > 0; i--) {
			*data = 0; data++;
		}
		data -= length * length;
	}
	matrix(int height, int width) {
		this->h = height, this->w = width;
		data = new T[h*w];
		for (int i = 0; i < h*w; i++) {
			*data = 0; data++;
		}
		data -= h * w;
	}
	matrix(int length, Matrix_type::_M_type type) {
		h = w = length;
		data = new T[length*length];
		for (int i = length * length; i > 0; i--) {
			*data = 0; data++;
		}
		data -= length * length;
		if (!type.type) return;
		if (type.type & 1) {
			for (int i = 0; i < length; i++) {
				*(data + i * w + i) = 1;
			}
			return;
		}
	}
	matrix(const matrix<T> &M) {
		h = M.h, w = M.w;
		if (h == 0 || w == 0) { data = 0; return; }
		data = new T[h * w];
		const T *p = M.data;
		for (int i = h * w; i > 0; i--) {
			*data = *p;
			data++, p++;
		}
		data -= h * w;
	}
	matrix(int h, int w, T p[]) {
		this->h = h, this->w = w;
		data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*data = *p; data++, p++;
		}
		data -= h * w, p -= h * w;
	}
	matrix(int length, T p[]) {
		this->h = this->w = length;
		data = new T[length*length];
		for (int i = length * length; i > 0; i--) {
			*data = *p; data++, p++;
		}
		data -= length * length, p -= length * length;
	}
	matrix(int h, int w, const initializer_list<T> &q) {
		this->h = h, this->w = w;
		data = new T[h*w];
		int ps = q.size();
		int ds = h * w;
		if (ps > ds) ps = ds;
		int i;
		const T* p = q.begin();
		for (i = 0; i < ps; i++) {
			*data = *p; data++, p++;
		}
		for (; i < ds; i++) {
			*data = *p; data++, p++;
		}
		data -= ds;
	}
	matrix(int length, const initializer_list<T> &q) {
		this->h = this->w = length;
		data = new T[h*w];
		int ps = q.size();
		int ds = h * w;
		if (ps > ds) ps = ds;
		int i;
		const T* p = q.begin();
		for (i = 0; i < ps; i++) {
			*data = *p; data++, p++;
		}
		for (; i < ds; i++) {
			*data = *p; data++, p++;
		}
		data -= ds;

	}
	matrix(const initializer_list<initializer_list<T>> &q) {
		h = q.size(); w = 0;
		if (h == 0) return;
		const initializer_list<T>* k = q.begin();
		while (k < q.end()) {
			if (k->size() > w) w = k->size();
			k++;
		}
		k -= h;
		const T* s;
		if (w == 0) return;
		data = new T[w*h];
		for (int i = 0; i < h; i++) {
			s = k->begin();
			for (int j = 0; j < w; j++) {
				if (s < k->end()) *data = *s;
				else *data = 0;
				data++, s++;
			}
			k++;
		}
		data -= h * w;
		return;
	}
	~matrix() {
		h = w = 0;
		delete data;
		data = 0;
	}
	matrix<T>& operator = (const matrix<T> &M) {
		const T* p = M.data;
		if (this->h != M.h || this->w != M.w) {
			this->h = M.h, this->w = M.w;
			if (data != 0) delete data;
			data = new T[h*w];
		}
		for (int i = h * w; i > 0; i--) {
			*data = *p;
			data++, p++;
		}
		data -= h * w;
		return *this;
	}
	matrix<T>& operator = (const initializer_list<initializer_list<T>> &q) {
		if (data != 0) delete data, data = 0;
		this->h = q.size(); w = 0;
		if (h == 0) { return *this; }
		const initializer_list<T>* k = q.begin();
		while (k < q.end()) {
			if (k->size() > w) w = k->size();
			k++;
		}
		k -= h;
		const T* s;
		if (w == 0) return *this;
		data = new T[w*h];
		for (int i = 0; i < h; i++) {
			s = k->begin();
			for (int j = 0; j < w; j++) {
				if (s < k->end()) *data = *s;
				else *data = 0;
				data++, s++;
			}
			k++;
		}
		data -= h * w;
		return *this;
	}
	inline T* operator [] (int k) const {
		return &data[k*w];
	}
	inline T& at(unsigned row, unsigned column) {
		if (row != 0 && column != 0) row--, column--;
		row %= h, column %= w;
		return *(data + row * w + column);
	}
	inline const unsigned height() const { return this->h; }
	inline const unsigned width() const { return this->w; }
	inline const unsigned size() const { return h * w; }
	inline const unsigned length() const { return h * w; }
	friend ostream& operator << (ostream& os, const matrix<T> &M)
	{
		/*os << "{ ";
		for (int i = 0; i < M.h; i++) {
			os << "{";
			for (int j = 0; j < M.w; j++) {
				os << M[i][j] << ",";
			}
			os << "\b}, ";
		}
		os << "\b\b }";
		return os;

		if (M.h == 0 || M.w == 0 || M.data == 0) {
			if (M.h == 0 && M.w == 0 && M.data == 0) {
				os << "#NAM" << endl;
			}
			else os << "<empty>" << endl;
			return os;
		}
		const T* p = M.data;
		for (int i = 0; i < M.h; i++) {
			for (int j = 0; j < M.w; j++) {
				os << *p << "\t";
				p++;
			}
			os << endl;
		}
		return os;*/
		os << "\\begin{bmatrix}";
		for (unsigned i = 0; i < M.height(); i++) {
			for (unsigned j = 0; j < M.width(); j++) {
				os << M[i][j] << "&";
			}
			os << "\b\\\\[3pt]";
		}
		os << "\\end{bmatrix}";
		return os;
	}

	inline bool empty() const {
		return h == 0 || w == 0 || data == 0;
	}
	inline bool square() const {
		return h == w && h != 0;
	}
	inline bool row() const {
		return h == 1 && w != 0;
	}
	inline bool column() const {
		return w == 1 && h != 0;
	}
	bool zero() {
		T* p = data;
		for (int i = h * w; i > 0; i--) {
			if (*p != 0) return 0;
		}
		return 1;
	}
	bool identity() {
		if (h != w || h == 0) return 0;
		T* p = data;
		T c = *p;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				if (i == j) {
					if (*p != 1) return 0;
				}
				else {
					if (*p != 0) return 0;
				}
				p++;
			}
		}
	}
	bool diagonal() {
		if (h != w || h == 0) return 0;
		T* p = data;
		T c = *p;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				if (i != j) {
					if (*p != 0) return 0;
				}
				else {
					if (*p != c) return 0;
				}
				p++;
			}
		}
		return 1;
	}
	bool upper_triangular() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < i; j++) {
				if ((*this)[i][j] != 0) return 0;
			}
		}
	}
	bool lower_triangular() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = i + 1; j < w; j++) {
				if ((*this)[i][j] != 0) return 0;
			}
		}
	}
	bool skew_upper_triangular() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = w - i; j < w; j++) {
				if ((*this)[i][j] != 0) return 0;
			}
		}
	}
	bool skew_lower_triangular() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w - i - 1; j++) {
				if ((*this)[i][j] != 0) return 0;
			}
		}
	}
	bool symmetric() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j <= i; j++) {
				if ((*this)[i][j] != (*this)[j][i]) return 0;
			}
		}
	}
	bool skew_symmetric() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j <= i; j++) {
				if ((*this)[i][j] != -(*this)[j][i]) return 0;
			}
			if ((*this)[i][i] != 0) return 0;
		}
	}
	bool invertible() {
		if (h != w || h == 0) return 0;
		if (h == 1) return *data != 0;
		if (h == 2) return (*data)*(*(data + 3)) - (*(data + 1)*(*(data + 2))) != 0;
		if (h == 3) return (*data)*((*(data + 4))*(*(data + 8)) - (*(data + 5)*(*(data + 7))))
			- (*(data + 1))*((*(data + 3))*(*(data + 8)) - (*(data + 5)*(*(data + 6))))
			+ (*(data + 2))*((*(data + 3))*(*(data + 7)) - (*(data + 4)*(*(data + 6)))) != 0;
		T *p, *q, c;
		matrix<T> A(*this);
		for (int i = 0; i < h; i++) {
			p = A.data + i * w + i;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) return 0;
					q += w;
					if (*q != 0) {
						for (int k = i; k < w; k++) {
							*p += *q; p++, q++;
						}
						p -= w - i;
						break;
					}
				}
			}
			q = p;
			for (int j = i + 1; j < h; j++) {
				q += w;
				if (*q != 0) {
					c = (*q) / (*p);
					for (int k = i; k < w; k++) {
						*q -= c * (*p);
						p++, q++;
					}
					p -= w - i, q -= w - i;
				}
			}
		}
		return 1;
	}
	bool elementary() {
		if (h != w) return 0;

	}
	inline bool col() { return column(); }
	inline bool sqr() { return square(); }
	inline bool utri() { return upper_triangular(); }
	inline bool ltri() { return lower_triangular(); }
	inline bool sutri() { return skew_upper_triangular(); }
	inline bool sltri() { return skew_lower_triangular(); }

	bool operator == (matrix<T> a) {
		if (h != a.h || w != a.w || a.h == 0) return 0;
		for (int i = h * w; i > 0; i--) {
			if (*data != *a.data) {
				data -= (h*w - i), a.data -= (a.h*a.w - i);
				return 0;
			}
			data++, a.data++;
		}
		data -= h * w, a.data -= a.h*a.w;
		return 1;
	}
	bool operator != (matrix<T> a) {
		if (h != a.h || w != a.w || a.h == 0 || a.w == 0) return 1;
		for (int i = h * w; i > 0; i--) {
			if (*data != *a.data) {
				data -= (h*w - i), a.data -= (a.h*a.w - i);
				return 1;
			}
			data++, a.data++;
		}
		data -= h * w, a.data -= a.h*a.w;
		return 0;
	}
	matrix<T> operator + (matrix<T> a) {
		if (h != a.h || w != a.w) {
			delete a.data; a.data = 0;
			a.h = a.w = 0;
			return a;
		}
		for (int i = h * w; i > 0; i--) {
			*a.data += *data; a.data++, data++;
		}
		data -= h * w, a.data -= h * w;
		return a;
	}
	void operator += (matrix<T> &a) {
		if (h != a.h || w != a.w) {
			delete data; data = 0;
			h = w = 0;
			return;
		}
		for (int i = h * w; i > 0; i--) {
			*data += *a.data; data++, a.data++;
		}
		data -= h * w, a.data -= h * w;
	}
	matrix<T> operator - () {
		matrix<T> a;
		a.h = h, a.w = w, a.data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*a.data = -*data;
			a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}
	matrix<T> operator - (matrix<T> a) {
		if (a.h != h || a.w != w) {
			delete a.data; a.data = 0;
			a.h = a.w = 0;
			return a;
		}
		for (int i = h * w; i > 0; i--) {
			*a.data = *data - *a.data; a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}
	void operator -= (matrix<T> &a) {
		if (a.h != h || a.w != w) {
			delete data; data = 0;
			h = w = 0;
			return;
		}
		for (int i = h * w; i > 0; i--) {
			*data -= *a.data; data++, a.data++;
		}
		data -= h * w, a.data -= h * w;
	}
	matrix<T> operator * (const T &b) {
		matrix<T> a;
		a.h = h, a.w = w, a.data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*a.data = *data * b;
			a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}
	void operator *= (const T &b) {
		for (int i = h * w; i > 0; i--) {
			*data *= b;
			data++;
		}
		data -= h * w;
	}
	matrix<T> operator / (const T &b) {
		matrix<T> a;
		a.h = h, a.w = w, a.data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*a.data = *data / b;
			a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}
	matrix<T> operator /= (const T &b) {
		for (int i = h * w; i > 0; i--) {
			*data /= b, data++;
		}
		data -= h * w;
	}
	/*template<typename t> matrix<T> operator * (const t &b) {
		matrix<T> a;
		a.h = h, a.w = w, a.data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*a.data = *data * b;
			a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}
	template<typename t> void operator *= (const t &b) {
		for (int i = h * w; i > 0; i--) {
			*data *= b;
			data++;
		}
		data -= h * w;
	}
	template<typename t> friend matrix<T> operator * (const t &b, matrix<T> a) {
		for (int i = a.h * a.w; i > 0; i--) {
			*a.data *= b;
			a.data++;
		}
		a.data -= a.h * a.w;
		return a;
	}
	template<typename t> matrix<T> operator / (const t &b) {
		matrix<T> a;
		a.h = h, a.w = w, a.data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*a.data = *data / b;
			a.data++, data++;
		}
		a.data -= h * w, data -= h * w;
		return a;
	}*/
	matrix<T> operator * (const matrix<T> &a) {
		/*if (h == 1 && w == 1) {
			for (int i = a.h*a.w; i > 0; i--) {
				*a.data *= *data; a.data++;
			}
			a.data -= a.h*a.w;
			return a;
		}
		if (a.h == 1 && a.w == 1) {
			for (int i = h * w; i > 0; i--) {
				*data *= *a.data, data++;
			}
			data -= h * w;
			return a;
		}*/
		if (w != a.h) {
			return matrix();
		}
		matrix b(h, w);
		for (int i = 0; i < b.h; i++) {
			for (int j = 0; j < b.w; j++) {
				for (int k = 0; k < w; k++) {
					b[i][j] += *(data + i * w + k) * (*(a.data + k * w + j));
				}
			}
		}
		return b;
	}
	void operator *= (matrix<T> a) {
		/*if (a.h == 1 && a.w == 1) {
			for (int i = h * w; i > 0; i--) {
				*data *= *a.data, data++;
			}
			data -= h * w;
			return;
		}*/
		if (w != a.h) {
			delete data; data = 0; h = w = 0; return;
		}
		matrix b(h, a.w);
		for (int i = 0; i < b.h; i++) {
			for (int j = 0; j < b.w; j++) {
				for (int k = 0; k < w; k++) {
					b[i][j] += *(data + i * w + k) * a[k][j];
				}
			}
		}
		delete data; h = b.h, w = b.w;
		data = new T[h*w];
		for (int i = h * w; i > 0; i--) {
			*data = *b.data; data++, b.data++;
		}
		data -= h * w, b.data -= b.h*b.w;
	}
	Vector<T> operator * (const Vector<T> &a) {
		if (a.l != this->w) return Vector<T>();
		Vector<T> b(this->h);
		for (int i = 0; i < h; i++) {
			b[i] = 0;
			for (int j = 0; j < w; j++) {
				b[i] += (*this)[i][j] * a[j];
			}
		}
		return b;
	}
	friend matrix<T> pow(matrix<T> a, int b) {
		if (a.h != a.w) {
			a.h = a.w = 0;
			delete a.data; a.data = 0;
			return a;
		}
		if (b == 0) {
			for (int i = a.h*a.w; i > 0; i--) {
				*a.data = 0; a.data++;
			}
			a.data -= a.h*a.w;
			for (int i = 0; i < a.h; i++) {
				a[i][i] = 1;
			}
			return a;
		}
		if (b < 0) {
			a = a.inverse();
			if (a.data == 0) return a;
			b = -b;
		}
		matrix<T> M(a.h, Matrix_type::_M_Identity), c(a.h);
		while (b != 0) {
			for (int i = 0; i < c.h; i++) {
				for (int j = 0; j < c.h; j++) {
					for (int k = 0; k < c.h; k++) {
						c[i][j] += M[i][k] * a[k][j];
					}
				}
			}
			for (int i = a.h*a.w; i > 0; i--) {
				*M.data = *c.data;
				*c.data = 0;
				M.data++, c.data++;
			}
			M.data -= c.h*c.h, c.data -= c.h*c.h;
			b--;
		}
		return M;
	}
	inline matrix<T> operator ^ (int b) {
		return pow(*this, b);
	}
	inline void operator ^= (int b) {
		*this = pow(*this, b);
	}
	matrix<T> transpose() {
		matrix<T> M(w, h);
		for (int i = 0; i < M.h; i++) {
			for (int j = 0; j < M.w; j++) {
				M[i][j] = (*this)[j][i];
			}
		}
		return M;
	}
	virtual matrix<T> inverse() {
		if (h != w || h == 0) return matrix();
		matrix<T> A(*this), I(h, Matrix_type::_M_Identity);
		long long dif = (long long)(I.data) - (long long)(A.data);
		T *p, *q, c;
		for (int i = 0; i < h; i++) {
			p = A.data + i * w + i;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) return matrix();
					q += w;
					if (*q != 0) {
						for (int k = i; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w, q -= w;
						*((unsigned*)&p) += dif, *((unsigned*)&q) += dif;
						for (int k = 0; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w - i, q -= w - i;
						*((unsigned*)&p) -= dif, *((unsigned*)&q) -= dif;
						break;
					}
				}
			}
			if (*p != 0) {
				if (*p != 1) {
					c = 1 / (*p);
					for (int j = i; j < w; j++) {
						*p *= c, p++;
					}
					p -= w, *((unsigned*)&p) += dif;
					for (int j = 0; j < w; j++) {
						*p *= c, p++;
					}
					p -= w - i, *((unsigned*)&p) -= dif;
				}
				for (int j = 0; j < h; j++) {
					p = A.data + i * w + i, q = A.data + j * w + i;
					if (j != i && *q != 0) {
						c = *q;
						for (int k = i; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
						p -= w, q -= w;
						*((unsigned*)&p) += dif, *((unsigned*)&q) += dif;
						for (int k = 0; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
					}
				}
			}
		}
		/*p = A.data;
		for (int i = 0; i < h; i++) {
			if (*p == 0) return matrix();
			p += w + 1;
		}*/
		return I;
	}
	virtual matrix<T> adjugate() {
		if (h != w || h == 0) return matrix();
		T *p, *q, c;
		matrix<T> M(*this), R(h, Matrix_type::_M_Identity);
		long long dif = (long long)(R.data) - (long long)(M.data);
		for (int i = 0; i < h; i++) {
			p = M.data + i * w + i;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) goto Not_Invertible;
					q += w;
					if (*q != 0) {
						for (int k = i; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w, q -= w;
						*((unsigned*)&p) += dif, *((unsigned*)&q) += dif;
						for (int k = 0; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w - i, q -= w - i;
						*((unsigned*)&p) -= dif, *((unsigned*)&q) -= dif;
						break;
					}
				}
			}
			if (*p != 0) {
				if (*p != 1) {
					c = 1 / (*p);
					for (int j = i; j < w; j++) {
						*p *= c, p++;
					}
					p -= w, *((unsigned*)&p) += dif;
					for (int j = 0; j < w; j++) {
						*p *= c, p++;
					}
					p -= w - i, *((unsigned*)&p) -= dif;
				}
				for (int j = 0; j < h; j++) {
					p = M.data + i * w + i, q = M.data + j * w + i;
					if (j != i && *q != 0) {
						c = *q;
						for (int k = i; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
						p -= w, q -= w;
						*((unsigned*)&p) += dif, *((unsigned*)&q) += dif;
						for (int k = 0; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
					}
				}
			}
		}
		c = this->Det();
		p = R.data;
		for (int i = h * w; i > 0; i--) {
			*p *= c, p++;
		}
		return R;
	Not_Invertible: {
		if (M.data != 0) delete M.data;
		M.h--, M.w--; M.data = new T[M.w*M.h];
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				p = data, q = M.data;
				for (int m = 0; m < i; m++) {
					for (int n = 0; n < w; n++) {
						if (n != j) *q = *p, q++;
						p++;
					}
				}
				p += w;
				for (int m = i + 1; m < h; m++) {
					for (int n = 0; n < w; n++) {
						if (n != j) *q = *p, q++;
						p++;
					}
				}
				*(R.data + j * w + i) = ((i + j) & 1) ? -M.Det() : M.Det();
			}
		}
		return R;
		}
	}
	inline void clear() {
		for (int i = h * w; i > 0; i--) {
			*data = 0, data++;
		}
		data -= h * w;
	}

	void interchange(unsigned i, unsigned j) {
		i--, j--;
		if (i >= h || j >= h) {
			delete data; h = w = 0; return;
		}
		if (i == j) return;
		T *p, *q, c;
		p = data + i * w, q = data + j * w;
		for (int n = 0; n < w; n++) {
			c = *p, *p = *q, *q = c;
			p++, q++;
		}
		p = q = 0;
	}
	void multiply(unsigned i, T n) {
		i--;
		if (i >= h) {
			delete data; h = w = 0; return;
		}
		T *p = data + i * w;
		for (int k = 0; k < w; k++) {
			*p *= n; p++;
		}
		p = 0;
	}
	void multiply(unsigned i, T n, unsigned j) {
		i--, j--;
		if (i >= h || j >= h) {
			delete data; h = w = 0; return;
		}
		T *p, *q;
		p = data + i * w, q = data + j * w;
		for (int k = 0; k < w; k++) {
			*q += n * (*p), q++, p++;
		}
		p = q = 0;
	}
	void col_interchange(unsigned i, unsigned j) {
		i--, j--;
		if (i >= w || j >= w) {
			delete data; h = w = 0; return;
		}
		if (i == j) return;
		T *p, *q, c;
		p = data + i, q = data + j;
		for (int n = 0; n < h; n++) {
			c = *p, *p = *q, *q = c;
			p += w, q += w;
		}
		p = q = 0;
	}
	void col_multiply(unsigned i, T n) {
		i--;
		if (i >= w) {
			delete data; h = w = 0; return;
		}
		T *p = data + i;
		for (int k = 0; k < h; k++) {
			*p *= n; p += w;
		}
		p = 0;
	}
	void col_multiply(unsigned i, T n, unsigned j) {
		i--, j--;
		if (i >= w || j >= w) {
			delete data; h = w = 0; return;
		}
		T *p, *q;
		p = data + i, q = data + j;
		for (int k = 0; k < h; k++) {
			*q += n * (*p);
			q += w, p += w;
		}
		p = q = 0;
	}
	inline void swaprow(unsigned i, unsigned j) { interchange(i, j); }
	inline void swapcol(unsigned i, unsigned j) { col_interchange(i, j); }
	virtual void elimination() {
		int y = 0;
		T *p, *q, c;
		for (int i = 0; i < h; i++) {
			//cout << "y: " << y << endl;
			p = data + i * w + y;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) {
						i--; break;
					}
					q += w;
					if (*q != 0) {
						for (int k = y; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w - y;
						//cout << "Interchange: " << endl << *this << endl << endl;
						break;
					}
				}
			}
			if (*p != 0) {
				if (*p != 1) {
					c = 1 / (*p);
					for (int j = y; j < w; j++) {
						*p *= c, p++;
					}
					//cout << "Multiply: " << endl << *this << endl << endl;
				}
				for (int j = 0; j < h; j++) {
					p = data + i * w + y, q = data + j * w + y;
					if (j != i && *q != 0) {
						c = *q;
						for (int k = y; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
						//cout << "Multiply + Add: " << endl << *this << endl << endl;
					}
				}
			}
			y++;
		}
	}
	virtual T Det() {
		if (h != w) return NAN;
		if (h == 0) return data == 0 ? 0 : 1;
		if (h == 1) return *data;
		if (h == 2) return (*data)*(*(data + 3)) - (*(data + 1)*(*(data + 2)));
		if (h == 3) return (*data)*((*(data + 4))*(*(data + 8)) - (*(data + 5)*(*(data + 7))))
			- (*(data + 1))*((*(data + 3))*(*(data + 8)) - (*(data + 5)*(*(data + 6))))
			+ (*(data + 2))*((*(data + 3))*(*(data + 7)) - (*(data + 4)*(*(data + 6))));
		matrix<T> A(*this);
		T c;
		/*T *p, *q;
		for (int i = 0; i < h; i++) {
			p = A.data + i * w + i;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) return 0;
					q += w;
					if (*q != 0) {
						for (int k = i; k < w; k++) {
							*p += *q; p++, q++;
						}
						p -= w - i;
						break;
					}
				}
			}
			q = p;
			for (int j = i + 1; j < h; j++) {
				q += w;
				if (*q != 0) {
					c = (*q) / (*p);
					for (int k = i; k < w; k++) {
						*q -= c * (*p);
						p++, q++;
					}
					p -= w - i, q -= w - i;
				}
			}
		}
		p = A.data; c = 1;
		for (int i = 0; i < h; i++) {
			if (*p == 0) return 0;
			c *= *p, p += w + 1;
		}*/
		for (int i = 0; i < h; i++) {
			for (int j = i + 1; j < h; j++) {
				c = -(A[j][i] / A[i][i]);
				for (int k = i; k < h; k++) A[j][k] += A[i][k] * c;
			}
		}
		c = 1;
		for (int i = 0; i < h; i++) c *= A[i][i];
		if (isnan(c)) return 0;
		return c;
	}
	friend inline T det(matrix &A) {
		return A.Det();
	}
	T Determinant() {
		if (h != w) return 0;
		if (h == 0) return data == 0 ? 0 : 1;
		if (h == 1) return *data;
		if (h == 2) return (*data)*(*(data + 3)) - (*(data + 1)*(*(data + 2)));
		if (h == 3) return (*data)*((*(data + 4))*(*(data + 8)) - (*(data + 5)*(*(data + 7))))
			- (*(data + 1))*((*(data + 3))*(*(data + 8)) - (*(data + 5)*(*(data + 6))))
			+ (*(data + 2))*((*(data + 3))*(*(data + 7)) - (*(data + 4)*(*(data + 6))));
		matrix a(h - 1);
		T *p, *q, c, det = 0;
		for (int i = 0; i < w; i++) {
			if (*(data + i) != 0) {
				p = data + w, q = a.data;
				for (int j = 1; j < h; j++) {
					for (int k = 0; k < w; k++) {
						if (k != i) *q = *p, q++;
						p++;
					}
				}
				if (i & 1) det -= *(data + i)*a.Determinant();
				else det += *(data + i)*a.Determinant();
			}
		}
		return det;
	}
	friend inline T determinant(matrix &A) {
		return A.Determinant();
	}
	friend inline matrix<T> adj(matrix &A) {
		return A.adjugate();
	}
	virtual unsigned rank() {
		//cout << *this << endl;
		matrix M(*this);
		int y = 0, a = w > h ? w : h;	// a: probably relating with rank
		bool it = (h == w);
		T *p, *q, c;
		for (int i = 0; i < h; i++) {
			//cout << "y: " << y << endl;
			p = M.data + i * w + y;
			if (*p == 0) {
				q = p;
				for (int j = i + 1; j <= h; j++) {
					if (j == h) {
						//cout << "Break once; \n";
						i--, a--; it = 0; break;
					}
					q += w;
					if (*q != 0) {
						for (int k = y; k < w; k++) {
							c = *q, *q = *p, *p = c;
							p++, q++;
						}
						p -= w - y;
						//cout << "Interchange: " << endl << M << endl << endl;
						break;
					}
				}
			}
			if (*p != 0) {
				for (int j = i + 1; j < h; j++) {
					p = M.data + i * w + y, q = M.data + j * w + y;
					if (j != i && *q != 0) {
						c = *q / (*p);
						for (int k = y; k < w; k++) {
							*q -= c * (*p);
							p++, q++;
						}
						//cout << "Multiply + Add: " << endl << M << endl << endl;
					}
				}
			}
			y++;
		}
		//cout << M << endl;
		if (it) return h;
		unsigned r = h > w ? w : h;
		for (r--; r >= 0; r--) {
			p = M.data + r * w;
			for (int i = 0; i < w; i++) {
				if (*p != 0) return r + 1;
				p++;
			}
		}
		return 0;
	}
	friend inline unsigned R(matrix &A) {
		return A.rank();
	}
	friend unsigned Rank(matrix A) {

	}

	matrix<T> submatrix(unsigned i, unsigned j) {
		i--, j--;
		if (i >= h || j >= w) return matrix();
		matrix<T> M(h - 1, w - 1);
		T *p = data, *q = M.data;
		for (int m = 0; m < i; m++) {
			for (int n = 0; n < w; n++) {
				if (n != j) *q = *p, q++;
				p++;
			}
		}
		p += w;
		for (int m = i + 1; m < h; m++) {
			for (int n = 0; n < w; n++) {
				if (n != j) *q = *p, q++;
				p++;
			}
		}
		return M;
	}
	matrix<T> submatrix(initializer_list<unsigned> rows, initializer_list<unsigned> cols) {
		vector<unsigned> r, c;
		for (const unsigned* i = rows.begin(); i < rows.end(); i++) {
			r.push_back((*i) - 1);
		}
		for (const unsigned* i = cols.begin(); i < cols.end(); i++) {
			c.push_back((*i) - 1);
		}
		bool sw = 0; int sz;
		sz = r.size();
		for (int i = 0; i < sz; i++) {
			if (r.at(i) > h) return matrix();
		}
		sz = c.size();
		for (int i = 0; i < sz; i++) {
			if (c.at(i) > h) return matrix();
		}
		do {
			sw = 0; sz = r.size();
			for (int i = 1; i < sz; i++) {
				if (r.at(i - 1) > r.at(i)) swap(r.at(i - 1), r.at(i)), sw = 1;
				else if (r.at(i - 1) == r.at(i)) r.erase(r.begin() + i), sz--, i--;
			}
		} while (sw);
		do {
			sw = 0; sz = c.size();
			for (int i = 1; i < sz; i++) {
				if (c.at(i - 1) > c.at(i)) swap(c.at(i - 1), c.at(i)), sw = 1;
				else if (c.at(i - 1) == c.at(i)) c.erase(c.begin() + i), sz--, i--;
			}
		} while (sw);
		if (r.size() > h || c.size() > w) return matrix();
		for (int i = 0; i < r.size(); i++) {
			cout << r.at(i) << " ";
		}
		cout << endl;
		for (int i = 0; i < c.size(); i++) {
			cout << c.at(i) << " ";
		}
		cout << endl;
		return matrix();
	}
	inline matrix<T> submat(unsigned i, unsigned j) { return submatrix(i, j); }

	inline T minor(unsigned i, unsigned j) {
		return submatrix(i, j).Det;
	}
	inline T cofactor(unsigned i, unsigned j) {
		return (i + j) & 1 ? -submatrix(i, j).Det : submatrix(i, j).Det;
	}
	// Minor, Cofactor...


	//friend void solve_q(matrix<T> A, Vector<T> &X, const Vector<T> &Y);
};

typedef matrix<double> Matrix;

// Solve systems of linear equations AX=Y quickly without checking. 
// Precondition: A must be a square matrix with non-zero determinant, X and Y must have been initialized. 
// Errors may occur when preconditions not satisfied. 
template<typename T> static void solve_q(matrix<T> A, Vector<T> &X, const Vector<T> &Y) {
	int l = Y.length();
	X = Y;
	T d;
	for (int i = 0; i < l; i++) {
		for (int j = i + 1; j < l; j++) {
			d = -(A[j][i] / A[i][i]);
			for (int k = i; k < l; k++) A[j][k] += A[i][k] * d;
			X[j] += X[i] * d;
		}
	}
	for (int i = l - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			d = -(A[j][i] / A[i][i]);
			for (int k = l - 1; k >= i; k--) A[j][k] += A[i][k] * d;
			X[j] += X[i] * d;
		}
	}
	for (int i = 0; i < l; i++) X[i] /= A[i][i];
}

#endif
