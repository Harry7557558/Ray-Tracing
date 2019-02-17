#include <iostream>
//#include <cmath>
#include <string>

using namespace std;

template <typename T> T abs(T a){
	if (a >= 0) return a;
	return -a;
}

class matrix;

class fraction{
	unsigned den;
	int num;
	static unsigned gcf(unsigned a, unsigned b)
	{
		unsigned c;
		if (a == 0 || b == 0) return 0;
		if (a == b) return a;
		if (b > a) c = a, a = b, b = c;
		c = 1;
		while (c != 0) c = a%b, a = b, b = c;
		return a;
	}
	static unsigned lcm(unsigned a, unsigned b)
	{
		unsigned long long c = a*b;
		if (c == 0) return 0;
		if (c == 1) return 1;
		return c / gcf(a, b);
	}
	static unsigned str2num(string a){
		unsigned n = 0;
		while (a.size() != 0){
			n *= 10, n += a[0] - 48;
			a.erase(0, 1);
		}
		return n;
	}
	static void rfcd(fraction &a, fraction &b)
	{
		unsigned c = lcm(a.den, b.den);
		a.num = c / a.den*a.num;
		b.num = c / b.den*b.num;
		a.den = b.den = c;
	}
	friend fraction operator ! (fraction a)
	{
		bool s = a.num < 0;
		fraction b;
		b.num = a.den, b.den = abs(a.num);
		if (s) b.num = -b.num;
		b.simplify();
		return b;
	}
public:
	fraction() :num(0), den(1) {/* cout << "constructor\n"; */}
	fraction(int a){
		num = a, den = 1;
	}
	fraction(double a){

	}
	fraction(string a){
		*this = a;
		return;

		bool sign = 0;
		if (a[0] == '-') sign = 1, a.erase(0, 1);
		if (a[a.find('/', 0) + 1] == '-') sign ^= 1;
		for (int i = 0; i < a.size(); i++){
			if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
		}
		if (a.find('/', 0) == -1){
			num = str2num(a), den = 1;
			if (sign) num = -num;
			return;
		}
		num = str2num(a.substr(0, a.find('/', 0)));
		den = str2num(a.substr(a.find('/', 0) + 1, a.size() - a.find('/', 0)));
		if (sign) num = -num;
		simplify();
	}
	inline friend ostream& operator << (ostream& os, const fraction &a_)
	{
		fraction a = a_;
		if (a.num == 0) {
			os << "0";
			return os;
		}
		if (a.den == 1){
			os << a.num;
			return os;
		}
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
			num = 0; return;
		}
		if (num == 0) {
			den = 1; return;
		}
		//cout << "Before:\t" << den << " " << num << endl;
		unsigned m = gcf(abs(num), den);
		if (m == 1) return;
		//cout << "After:\t" << den << " " << num << " " << m << endl;
		bool s = num < 0;
		den /= m; num = abs(num) / m;
		if (s) num = -num;
		//cout << "Finary:\t" << den << " " << num << " " << m << endl;
	}
#ifdef Debug_Succeed
	friend fraction operator / (const int& a, const int& b){
		fraction s;
		s.num = abs(a), s.den = abs(b);
		if ((a < 0) ^ (b < 0)) s.num = -s.num;
		s.simplify();
		return s;
	}
#endif
	fraction& operator = (int a){
		num = a, den = 1;
		return *this;
	}
	fraction& operator = (double a){
		fraction f;
		int w = a;
		if (a == w){
			f.num = w, f.den = 1;
			return f;
		}
		unsigned n = 1, d = 2;
		double s;
		do{
			s = (w*d + n) / d;
			if (s > a){
				d *= 2, n = n * 2 - 1;
				
			}
			if (s < a){
				d *= 2, n = n * 2 + 1;
			}
		} while (s != a);
	}
	void operator = (string a){
		bool sign = 0;
		if (a[0] == '-') sign = 1, a.erase(0, 1);
		if (a[a.find('/', 0) + 1] == '-' && a.find('/', 0) != 0) sign ^= 1;
		if (a.find('/', 0) == -1){
			if (a.find('.', 0) != -1){
				for (int i = 0; i < a.size(); i++){
					if ((a[i] < 47 || a[i] > 57) && a[i] != '.') a.erase(i, 1);
				}
				den = 1;
				for (int i = a.size()-a.find('.', 0)-1; i > 0; i--){
					den *= 10;
				}
				a.erase(a.find('.', 0), 1);
				num = str2num(a);
				if (sign) num = -num;
				simplify();
				return;
			}
			for (int i = 0; i < a.size(); i++){
				if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
			}
			num = str2num(a), den = 1;
			if (sign) num = -num;
			simplify();
			return;
		}
		//cout << a << endl;
		for (int i = 0; i < a.size(); i++){
			if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
		}
		num = str2num(a.substr(0, a.find('/', 0)));
		den = str2num(a.substr(a.find('/', 0) + 1, a.size() - a.find('/', 0)));
		if (sign) num = -num;
		simplify();
	}
	inline friend bool operator == (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		a.simplify(), b.simplify();
		if (a.den == b.den && a.num == b.num) return 1;
		return 0;
	}
	inline friend bool operator != (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		if (a == b) return 0;
		return 1;
	}
	inline friend bool operator == (fraction a, int b)
	{
		if (a.den == 0) return 0;
		if (b*a.den == a.num) return 1;
		return 0;
	}
	inline friend bool operator != (fraction a, int b)
	{
		if (a.den == 0) return 0;
		if (b*a.den == a.num) return 0;
		return 1;
	}
	inline friend bool operator < (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		rfcd(a, b);
		if (a.num < b.num) return 1;
		return 0;
	}
	inline friend bool operator >(fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0) return 0;
		rfcd(a, b);
		if (a.num > b.num) return 1;
		return 0;
	}
	inline friend bool operator <= (fraction a, fraction b)
	{
		if (a == b || a < b) return 1;
		return 0;
	}
	inline friend bool operator >= (fraction a, fraction b)
	{
		if (a == b || a > b) return 1;
		return 0;
	}
	inline friend fraction operator + (fraction a, int b)
	{
		a.num += b*a.den;
		a.simplify();
		return a;
	}
	inline friend void operator += (fraction &a, int b)
	{
		a.num += b*a.den;
		a.simplify();
	}
	inline friend fraction operator + (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0){
			a.den = a.num = 0; return a;
		}
		if (a.num == 0) {
			b.simplify(); return b;
		}
		if (b.num == 0){
			a.simplify(); return a;
		}
		rfcd(a, b);
		a.num += b.num;
		a.simplify();
		return a;
	}
	inline friend void operator += (fraction &a, fraction b)
	{
		if (a.den == 0 || b.den == 0){
			a.den = a.num = 0; return;
		}
		if (b.num == 0){
			a.simplify(); return;
		}
		if (a.num == 0) {
			b.simplify(); a = b; return;
		}
		rfcd(a, b);
		a.num += b.num;
		a.simplify();
	}
	inline friend fraction operator - (fraction a)
	{
		a.simplify();
		a.num = -a.num;
		return a;
	}
	inline friend fraction operator - (fraction a, int b)
	{
		return a + (-b);
	}
	inline friend void operator -= (fraction &a, int b)
	{
		a += (-b);
	}
	inline friend fraction operator - (fraction a, fraction b)
	{
		if (a.den == 0 || b.den == 0){
			a.den = a.num = 0; return a;
		}
		if (a.num == 0) {
			return -b;
		}
		if (b.num == 0){
			a.simplify(); return a;
		}
		rfcd(a, b);
		a.num -= b.num;
		a.simplify();
		return a;
	}
	inline friend void operator -= (fraction &a, fraction b)
	{
		if (a.den == 0 || b.den == 0){
			a.den = a.num = 0; return;
		}
		if (b.num == 0){
			a.simplify(); return;
		}
		if (a.num == 0) {
			b.simplify(); a = -b; return;
		}
		rfcd(a, b);
		a.num -= b.num;
		a.simplify();
	}
	inline friend fraction operator * (fraction a, int b)
	{
		a.num *= b;
		a.simplify();
		return a;
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
	inline friend fraction operator ^ (fraction a, unsigned k)
	{
		fraction b; b.den = 1, b.num = 1;
		a.simplify();
		for (int i = 0; i < k; i++){
			b.den *= a.den;
			b.num *= a.num;
		}
		return b;
	}
	inline friend void operator ^= (fraction &a, unsigned k)
	{
		a = a^k;
	}

	friend class matrix;
	friend ostream& operator << (ostream& os, const matrix &a);
	friend bool operator == (const matrix &a, const matrix &b);
	friend bool operator != (const matrix &a, const matrix &b);
	friend matrix operator + (const matrix &a, const matrix &b_);
	friend void operator += (matrix &a, const matrix &b);
	friend matrix operator - (const matrix &a);
	friend matrix operator - (const matrix &a, const matrix &b_);
	friend void operator -= (matrix &a, const matrix &b);
	friend matrix operator * (matrix a, fraction k);
	friend void operator *= (matrix &a, fraction k);
	friend matrix operator * (matrix a, int k);
	friend void operator *= (matrix &a, int k);
	friend matrix operator * (const matrix &a_, const matrix &b_);
	friend void operator *= (matrix &a, const matrix &b);
	friend matrix operator ^ (const matrix &a_, int k);
	friend void operator ^= (matrix &a, int k);
	friend void operator >>= (matrix &a, int k);
	friend void operator <<= (matrix &a, int k);
	friend matrix operator >> (const matrix &a_, int k);
	friend matrix operator << (const matrix &a_, int k);
	friend class polynomial;
};

#ifdef It_Still_Use_in_Class_Fraction
inline bool OutputFraction(fraction a_)
{
	fraction a = a_;
	if (a.den == 0){
		cout << "ERROR!\n";
		return 0;
	}
	if (a.den < 0) a.den = -a.den, a.num = -a.num;
	cout << a.num << "/" << a.den;
	return 1;
}
bool simple(fraction &a)
{
	if (a.den < 0) a.den = -a.den, a.num = -a.num;
	bool s = 0;
	if (a.num < 0) s = 1, a.num = -a.num;
	if (a.den == 0) {
		a.num = 0; return 0;
	}
	if (a.num == 0) {
		a.den = 1; return 1;
	}
	int m = gcd(a.num, a.den);
	a.den /= m, a.num /= m;
	if (s == 1) a.num = -a.num;
	return 1;
}
#endif

