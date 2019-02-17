/*
Large Integer
Largest: 2¹⁶³⁸⁴-1   1,189,731,495,357,231,765,085,759,326,628,007... nearing 5000 bits. Almost infinity. 
*/

#include <iostream>
#include <cmath>
using namespace std; 
typedef unsigned long long gint; 
typedef unsigned char byte; 

/*
length: 16
used: 5
0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
xx xx xx xx xx 00 00 00 00 00 00 00 00 00 00 00
number = [0] + [1]×2⁶⁴ + [2]×2¹²⁸ + [3]×2²⁵⁶ + [4]×2⁵¹²
*/

class integer
{
	unsigned length, used; 
	bool sign; 
	gint *data; 
	void setwidth(unsigned n){
		gint *p = data;
		data = new gint[n];
		for (int i = 0; i < used; i++){
			*data = *p;
			data++, p++;
		}
		p -= used; delete p;
		for (int i = n - used; i > 0; i--){
			*data = 0, data++;
		}
		data -= n;
		length = n;
	}
	friend integer operator & (integer &a, integer &b){
		cout << "========= And Begin =========" << endl;
		integer c;
		cout << c.data << endl;
		cout << "Ping1\n";
		c.used = a.used > b.used ? b.used : a.used;
		if (c.used > c.length) c.setwidth(a.length > b.length ? b.length : a.length);
		c.data += c.used, b.data += c.used, a.data += c.used;
		cout << "Ping2\n";
		for (int i = 0; i < c.used; i++){
			c.data--, b.data--, a.data--;
			*c.data = (*a.data) & (*b.data);
		}
		cout << "Ping3\n";
		cout << "========== And End ==========" << endl;
		return c;
	}
	friend integer operator | (integer &a, integer &b){
		cout << "========= And Begin =========" << endl;
		integer c;
		cout << c.data << endl;
		cout << "Ping1\n";
		c.used = a.used < b.used ? b.used : a.used;
		if (c.used > c.length) c.setwidth(a.length < b.length ? b.length : a.length);
		c.data += c.used, b.data += c.used, a.data += c.used;
		cout << "Ping2\n";
		for (int i = 0; i < c.used; i++){
			c.data--, b.data--, a.data--;
			*c.data = (*a.data) | (*b.data);
		}
		cout << "Ping3\n";
		cout << "========== And End ==========" << endl;
		return c;
	}
	friend integer operator ^ (integer &a, integer &b){
		cout << "========= And Begin =========" << endl;
		integer c;
		cout << c.data << endl;
		cout << "Ping1\n";
		c.used = a.used < b.used ? b.used : a.used;
		if (c.used > c.length) c.setwidth(a.length < b.length ? b.length : a.length);
		c.data += c.used, b.data += c.used, a.data += c.used;
		cout << "Ping2\n";
		for (int i = 0; i < c.used; i++){
			c.data--, b.data--, a.data--;
			*c.data = (*a.data) ^ (*b.data);
		}
		cout << "Ping3\n";
		cout << "========== And End ==========" << endl;
		return c;
	}
public: 
	integer(){
		cout << "===== Constructor Begin =====" << endl;
		length = 4, used = 1, sign = 0;
		data = new gint[4];
		data[0] = data[1] = data[2] = data[3] = 0;
		cout << "====== Constructor End ======" << endl;
	}
	template <typename T> integer(T a){
		cout << "==== CpConstructor Begin ====" << endl;
		length = 4, used = 1, sign = 0;
		data = new gint[4];
		data[0] = data[1] = data[2] = data[3] = 0;
		if (a < 0) sign = 1, a = -a;
		else sign = 0;
		for (int i = 0; i < used; i++){
			*data = 0, data++;
		}
		data -= used;
		*data = a;
		if (sign) a = -a;
		cout << "===== CpConstructor End =====" << endl;
	}
	integer(integer &a){
		cout << "==== CPConstructor Begin ====" << endl;
		length = 4, used = 1, sign = 0;
		data = new gint[4];
		data[0] = data[1] = data[2] = data[3] = 0;
		if (length < a.used) setwidth(a.length);
		used = a.used;
		for (int i = 0; i < used; i++){
			*data = *a.data;
			data++, a.data++;
		}
		data -= used, a.data -= used;
		cout << "===== CPConstructor End =====" << endl;
	}
	~integer(){
		cout << "===== Destructor Begin =====" << endl;
		cout << data << endl;
		if (data != 0) delete data;
		cout << data << endl;
		data = 0;
		cout << "====== Destructor End ======" << endl;
	}
	template <typename T> void operator = (T &a){
		cout << "======= CpyInt Begin =======" << endl;
		if (a < 0) sign = 1, a = -a;
		else sign = 0;
		data += used;
		for (int i = 0; i < used; i++){
			data--, *data = 0;
		}
		*data = a;
		used = 1;
		cout << "======== CpyInt End ========" << endl;
	}
	void operator = (integer &a){
		cout << "======= Copying Begin ======" << endl;
		sign = a.sign;
		if (length < a.used) setwidth(a.length);
		used = a.used;
		for (int i = 0; i < used; i++){
			*data = *a.data;
			data++, a.data++;
		}
		data -= used, a.data -= used;
		cout << "======= Copying End. =======" << endl;
	}
	friend istream& operator >> (istream& is, integer &a){
		string n; cin >> n;
		if (n[0] == '-') a.sign = 1; else a.sign = 0;
		for (int i = n.size() - 1; i >= 0; i--){
			if (n[i] < '0' || n[i] > '9') n.erase(i, 1);
		}
	}
	friend ostream& operator << (ostream& os, integer &a){
		cout << "======= Output Begin =======" << endl;
		if (a.sign) os << "-";
		a.data += a.used - 1;
		for (int i = a.used; i > 0; i--){
			for (int j = 63; j >= 0; j--){
				os << (*a.data >> j & 1);
			}
			a.data--;
		}
		a.data++;
		cout << "\n======== Output End ========" << endl;
		return os;
	}
};

#ifdef USED
class integer
{
	gint *data; 
	inline friend integer operator & (const integer &a, const integer &b){
		integer c; 
		gint *p = c.data, *s = a.data, *t = b.data;
		for (int i = 0; i < 256; i++){
			*p = *s & *t;
			p++, s++, t++; 
		}
		return c;
	}
	inline friend void operator &= (integer &a, const integer &b){
		gint *p = a.data, *q = b.data; 
		for (int i = 0; i < 256; i++){
			*p &= *q;
			p++, q++; 
		}
	}
	inline friend integer operator | (const integer &a, const integer &b){
		integer c;
		gint *p = c.data, *s = a.data, *t = b.data;
		for (int i = 0; i < 256; i++){
			*p = *s | *t;
			p++, s++, t++;
		}
		return c;
	}
	inline friend void operator |= (integer &a, const integer &b){
		gint *p = a.data, *q = b.data;
		for (int i = 0; i < 256; i++){
			*p |= *q;
			p++, q++;
		}
	}
	inline friend integer operator ^ (const integer &a, const integer &b){
		integer c;
		gint *p = c.data, *s = a.data, *t = b.data;
		for (int i = 0; i < 256; i++){
			*p = *s ^ *t;
			p++, s++, t++;
		}
		return c;
	}
	inline friend void operator ^= (integer &a, const integer &b){
		gint *p = a.data, *q = b.data;
		for (int i = 0; i < 256; i++){
			*p ^= *q;
			p++, q++;
		}
	}
	inline friend integer operator ~ (const integer &a){
		integer b; 
		gint *p = a.data, *q = b.data;
		for (int i = 0; i < 256; i++){
			*p = ~*q;
			p++, q++;
		}
		return b;
	}
	inline friend integer operator << (const integer &a, int k){
		if (!k) return a;
		if (k < 0) return a >> -k;
		integer b;
		gint *p = a.data + 256, *q = b.data + 256;
		for (int i = 0; i < 256; i++){
			q--;
			*q = *p >> (64 - k);
			p--;
			*q |= *p << k;
		}
		p = b.data + 255;
		*p = *p >> k << k;
		return b;
	}
	inline friend void operator <<= (integer &a, int k){
		a = a << k;
	}
	inline friend integer operator >> (const integer &a, int k){
		if (!k) return a;
		if (k < 0) return a << -k;
		integer b;
		gint *p = a.data + 255, *q = b.data + 256;
		for (int i = 0; i < 256; i++){
			q--;
			*q = *p >> k;
			p--;
			*q |= *p << (64 - k);
		}
		*q = *q << k >> k;
		return b;
	}
	inline friend void operator >>= (integer &a, int k){
		a = a >> k;
	}
	// In << and >> operation, k should less than 64; 
public:
	integer(){
		data = new gint[256];
		gint *p = data;
		for (int i = 0; i < 256; i++){
			*p = 0;
			p++;
		}
	}
	integer(const integer &a){
		gint *p = data, *q = a.data;
		for (int i = 0; i < 256; i++){
			*p = *q;
			p++, q++;
		}
	}
	~integer(){
		delete data;
		data = 0;
	}
	integer& operator = (unsigned long long &a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (unsigned a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (unsigned short a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (unsigned char a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (long long a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this; 
	}
	integer& operator = (int a){
		bool s = a >> 32;
		if (s) a = -a;
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (short a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (char a){
		gint *p = data;
		for (int i = 0; i < 256; i++) {
			*p = 0, p++;
		}
		*(p - 1) = a;
		return *this;
	}
	integer& operator = (integer a){
		gint *p = data, *q = a.data;
		for (int i = 0; i < 256; i++){
			*p = *q;
			p++, q++;
		}
		return *this; 
	}
	friend integer operator + (const integer &a_, const integer &b_){
		integer a = a_, b = b_, t;
		/*while (b){
			t = a, a = t^b, b = (t&b) << 1;
		}*/
		return a;
	}
#include <iomanip>
	friend ostream& operator << (ostream& os, const integer &a){
		cout << &a << endl;
		byte const *p = (byte*)&(*a.data);
		cout << p << endl;
		os << ((*p) >> 7 & 1) << endl;
		for (int i = 0; i < 1024; i++){
			for (int j = 7; j >= 0; j--){
				os << ((*p) >> j & 1);
			}
			p++;
		}
		os << endl;
		return os;
	}
}; 
#endif