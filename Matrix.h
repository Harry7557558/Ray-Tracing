﻿#pragma once
#pragma warning(disable: 4244)

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <queue>
#include <stack>

#include <cstdlib>
#include <ctime>

using namespace std;

#ifndef LargeInteger

//#define LargeInteger

class integer {
	vector<unsigned> data;
	bool sign;
	// take care not to place an empty object at the back of data. 
	// if the value of a integer is zero, the sign should be zero, and there's only one empty value in data domain. 
public:
	integer() : sign(0) {
		data.clear();
		data.push_back(0);
	}
	integer(const int &a) {
		data.clear();
		data.push_back(abs(a));
		sign = a < 0;
	}
	integer(const unsigned &a) {
		data.clear();
		data.push_back(a);
		sign = 0;
	}
	integer(long long a) {
		if (a < 0) sign = 1, a = -a;
		else sign = 0;
		data.clear();
		data.push_back(a);
		a >>= 32;
		if (!a) data.push_back(a);
	}
	integer(unsigned long long a) {
		data.clear();
		data.push_back(a);
		a >>= 32;
		if (!a) data.push_back(a >> 32);
	}
	integer(const double& fb) {
		unsigned long long a = *((unsigned long long*)&fb);
		sign = a >> 63;
		a <<= 1, a >>= 1;
		data.clear();
		unsigned exp = a >> 52;
		if (exp < 1023) {
			sign = 0;
			data.push_back(0);
			return;
		}
		exp -= 1023;
		if (!exp) {
			data.push_back(1);
			return;
		}
		if (exp < 64) {
			unsigned long long k = fb;
			data.push_back(k);
			data.push_back(k >> 32);
			return;
		}
		a <<= 12, a >>= 12;
		if (exp <= 52) {
			a >>= (52 - exp);
			a |= ((unsigned long long)1) << exp;
			data.push_back(a);
			a >>= 32;
			if (a) data.push_back(a);
			return;
		}

		for (int i = 63; i >= 0; i--) cout << ((a&(((unsigned long long)1) << i)) >> i);
		cout << endl;

		data.resize(exp >> 5, 0);
		data.push_back((1 << (exp & 31)) | (a >> (52 - (exp & 31))));
		if (exp & 31 >= 20) {
			data.at(data.size() - 2) = a << ((exp & 31) - 20);
			return;
		}
		data.at(data.size() - 2) = a >> (20 - (exp & 31));

		for (int i = 31; i >= 0; i--) cout << (((data.at(data.size() - 2))&(((unsigned)1) << i)) >> i);
		cout << endl;

		data.at(data.size() - 3) = a << ((exp & 31) + 12);

		for (int i = 31; i >= 0; i--) cout << (((data.at(data.size() - 3))&(((unsigned)1) << i)) >> i);
		cout << endl;

		/*
		A debug part in function main:

		integer a; double b;
		unsigned long long k;
		while (1) {
			cin >> b;
			k = *((unsigned long long*)&b);
			cout << ((k&(((unsigned long long)1) << 63)) >> 63) << " ";
			for (int i = 62; i >= 52; i--) cout << ((k&(((unsigned long long)1) << i)) >> i);
			cout << " ";
			for (int i = 51; i >= 0; i--) cout << ((k&(((unsigned long long)1) << i)) >> i);
			cout << endl;
			a = b;
			cout << a << endl << endl;
		}


		*/
	}
	integer(const integer& a) {
		data.resize(a.data.size());
		sign = a.sign;
		for (int i = 0; i < a.data.size(); i++) {
			data.at(i) = a.data.at(i);
		}
	}
	integer(string s) {
		if (s[0] == '-') sign = 1;
		else sign = 0;
		data.clear();
		if (s.find('.', 0) != -1) {
			s.erase(s.find('.', 0), s.size() - s.find('.', 0));
		}
		for (int i = 0; i < s.size(); i++) {
			if (s[i] < '0' || s[i] > '9') s.erase(i, 1), i--;
		}
		while (s[0] == '0') s.erase(0, 1);
		if (s.empty()) {
			data.push_back(0);
			sign = 0;
			return;
		}
		queue<bool> b;
		bool sg;
		string t;
		while (!s.empty()) {
			if (s[s.size() - 1] & 1) {
				s[s.size() - 1]--; b.push(1);
			}
			else b.push(0);
			sg = 0;
			while (!s.empty()) {
				if (sg) s[0] += 10;
				sg = s[0] & 1;
				t += "0"; t[t.size() - 1] = (s[0] + '0') / 2;
				s.erase(0, 1);
			}
			s = t, t.clear();
			while (s[0] == '0') s.erase(0, 1);
		}
		unsigned c;
		stack<bool> sb;
		while (b.size()) {
			for (int i = 0; i < 32; i++) {
				sb.push(b.front()), b.pop();
				if (!b.size()) break;
			}
			c = 0;
			while (sb.size()) {
				c <<= 1, c |= sb.top();
				sb.pop();
			}
			data.push_back(c);
			if (b.size() == 0) break;
		}
	}
	inline integer& operator = (const int &a) {
		data.clear();
		data.push_back(abs(a));
		sign = a < 0;
		return *this;
	}
	inline integer& operator = (const unsigned &a) {
		data.clear();
		data.push_back(a);
		sign = 0;
		return *this;
	}
	inline integer& operator = (long long a) {
		if (a < 0) sign = 1, a = -a;
		else sign = 0;
		data.clear();
		data.push_back(a);
		a >>= 32;
		if (!a) data.push_back(a);
		return *this;
	}
	inline integer& operator = (unsigned long long a) {
		data.clear();
		data.push_back(a);
		a >>= 32;
		if (!a) data.push_back(a >> 32);
		return *this;
	}
	inline integer& operator = (const integer& a) {
		data.resize(a.data.size());
		sign = a.sign;
		for (int i = 0; i < a.data.size(); i++) {
			data.at(i) = a.data.at(i);
		}
		return *this;
	}

	friend istream& operator >> (istream& is, integer &a) {
		string s;
		is >> s;

		//clock_t time = clock();

		if (s[0] == '-') a.sign = 1;
		else a.sign = 0;
		a.data.clear();
		if (s.find('.', 0) != -1) {
			s.erase(s.find('.', 0), s.size() - s.find('.', 0));
		}
		for (int i = 0; i < s.size(); i++) {
			if (s[i] < '0' || s[i] > '9') s.erase(i, 1), i--;
		}
		while (s[0] == '0') s.erase(0, 1);
		if (s.empty()) {
			a.data.push_back(0);
			a.sign = 0;
			return is;
		}

		queue<bool> b;
		bool sg;
		string t;
		while (!s.empty()) {
			if (s[s.size() - 1] & 1) {
				s[s.size() - 1]--; b.push(1);
			}
			else b.push(0);
			sg = 0;
			while (!s.empty()) {
				if (sg) s[0] += 10;
				sg = s[0] & 1;
				t += "0"; t[t.size() - 1] = (s[0] + '0') / 2;
				s.erase(0, 1);
			}
			s = t, t.clear();
			while (s[0] == '0') s.erase(0, 1);
		}

		unsigned c;
		stack<bool> sb;
		while (b.size()) {
			for (int i = 0; i < 32; i++) {
				sb.push(b.front()), b.pop();
				if (!b.size()) break;
			}
			c = 0;
			while (sb.size()) {
				c <<= 1, c |= sb.top();
				sb.pop();
			}
			a.data.push_back(c);
			if (b.size() == 0) break;
		}

		//cout << "Used " << (clock() - time) << "ms to Input.\n";

		return is;
	}
	friend ostream& operator << (ostream& os, const integer &a) {
		if (a.data.empty()) {
			os << "0";
			return os;
		}
		if (a.data.size() == 1) {
			if (!a.data.at(0)) {
				os << "0"; return os;
			}
			if (a.sign) os << "-";
			os << a.data.at(0);
			return os;
		}
		if (a.sign) os << "-";
		if (a.data.size() == 2) {
			unsigned long long g = a.data.at(1);
			g <<= 32; g |= a.data.at(0);
			os << g;
			return os;
		}

		//clock_t time = clock();

		stack<bool> b;
		unsigned u;
		for (int i = 0; i < a.data.size(); i++) {
			u = a.data.at(i);
			for (int j = 0; j < 32; j++) {
				b.push(u & 1), u >>= 1;
			}
		}
		while (!b.top()) b.pop();
		string s = "0";
		int n;
		bool sg;
		while (!b.empty()) {
			for (n = s.size() - 1; n >= 0; n--) {
				s[n] = s[n] * 2 - '0';
				if (s[n + 1] > '9') {
					s[n + 1] -= 10;
					s[n]++;
				}
			}
			if (s[0] > '9') {
				s[0] -= 10;
				s.insert(0, "1", 0, 1);
			}
			if (b.top()) {
				n = s.size() - 1;
				s[n]++;
				while (s[n] > '9') {
					s[n] -= 10, n--;
					if (n == -1) {
						s.insert(0, "1", 0, 1);
						break;
					}
					s[n]++;
				}
			}
			b.pop();
		}
		while (s[0] == '0') s.erase(0, 1);
		if (s.empty()) s = "0";
		os << s;

		//cout << "Used " << (clock() - time) << "ms to Output.\n";

		return os;
	}

	inline bool empty() {
		// F1. See if the value of the integer is equal to 0; 
		// F2. Format the integer; 
		if (data.empty()) {
			data.push_back(0); sign = 0;
			return 1;
		}
		while (!data.back()) {
			data.pop_back();
			if (data.empty()) {
				data.push_back(0); sign = 0;
				return 1;
			}
		}
		return 0;
	}
	inline bool negative() {
		return sign;
	}
	inline bool positive() {
		return sign && (data.size() == 1) && (!data.front());
	}
	inline bool odd() {
		return data.front() & 1;
	}
	inline bool even() {
		return !(data.front() & 1);
	}

	friend inline integer abs(integer a) {
		a.sign = 0; return a;
	}
	inline void Abs() {
		sign = 0;
	}
	inline void minus() {
		sign ^= 1;
	}
	inline void clear() {
		data.clear(); data.push_back(0); sign = 0;
	}
	//friend inline void swap(integer &a, integer &b);

	inline bool operator == (const integer &a) {
		if (data.size() != a.data.size()) return 0;
		if (sign != a.sign) return 0;
		for (int i = 0; i < data.size(); i++) {
			if (data.at(i) != a.data.at(i)) return 0;
		}
		return 1;
	}
	inline bool operator != (const integer &a) {
		if (data.size() != a.data.size()) return 1;
		if (sign != a.sign) return 1;
		for (int i = 0; i < data.size(); i++) {
			if (data.at(i) != a.data.at(i)) return 1;
		}
		return 0;
	}
	inline bool operator > (const integer &a) {
		if (sign != a.sign) return a.sign;
		if (sign) {
			if (data.size() != a.data.size()) return data.size() < a.data.size();
			return data.back() < a.data.back();
		}
		if (data.size() != a.data.size()) return data.size() > a.data.size();
		return data.back() < a.data.back();
	}
	inline bool operator >= (const integer &a) {
		if (sign != a.sign) return a.sign;
		if (sign) {
			if (data.size() != a.data.size()) return data.size() < a.data.size();
			for (int i = data.size() - 1; i >= 0; i--) {
				if (data.at(i) > a.data.at(i)) return 0;
				if (data.at(i) < a.data.at(i)) return 1;
			}
			return 1;
		}
		if (data.size() != a.data.size()) return data.size() > a.data.size();
		for (int i = data.size() - 1; i >= 0; i--) {
			if (data.at(i) > a.data.at(i)) return 1;
			if (data.at(i) < a.data.at(i)) return 0;
		}
		return 1;
	}
	inline bool operator < (const integer &a) {
		if (sign != a.sign) return sign;
		if (sign) {
			if (data.size() != a.data.size()) return data.size() > a.data.size();
			return data.back() > a.data.back();
		}
		if (data.size() != a.data.size()) return data.size() < a.data.size();
		return data.back() > a.data.back();
	}
	inline bool operator <= (const integer &a) {
		if (sign != a.sign) return sign;
		if (sign) {
			if (data.size() != a.data.size()) return data.size() > a.data.size();
			for (int i = data.size() - 1; i >= 0; i--) {
				if (data.at(i) < a.data.at(i)) return 0;
				if (data.at(i) > a.data.at(i)) return 1;
			}
			return 1;
		}
		if (data.size() != a.data.size()) return data.size() < a.data.size();
		for (int i = data.size() - 1; i >= 0; i--) {
			if (data.at(i) < a.data.at(i)) return 1;
			if (data.at(i) > a.data.at(i)) return 0;
		}
		return 1;
	}

	inline bool operator == (const int &a) {
		if (sign != (a < 0)) return 0;
		if (data.size() > 1) return 0;
		return data.front() == abs(a);
	}
	inline bool operator == (const unsigned &a) {
		if (sign) return 0;
		if (data.size() > 1) return 0;
		return data.front() == a;
	}
	inline bool operator != (const int &a) {
		if (sign != (a < 0)) return 1;
		if (data.size() > 1) return 1;
		return data.front() != abs(a);
	}
	inline bool operator != (const unsigned &a) {
		if (sign) return 1;
		if (data.size() > 1) return 1;
		return data.front() != a;
	}
	inline bool operator > (const int &a) {
		if (sign != (a < 0)) return a < 0;
		if (sign) {
			if (data.size() > 1) return 0;
			return data.front() < abs(a);
		}
		if (data.size() > 1) return 1;
		return data.front() > abs(a);
	}
	inline bool operator > (const unsigned &a) {
		if (sign) return 0;
		if (data.size() > 1) return 1;
		return data.front() > a;
	}
	inline bool operator >= (const int &a) {
		if (sign != (a < 0)) return a < 0;
		if (sign) {
			if (data.size() > 1) return 0;
			return data.front() <= abs(a);
		}
		if (data.size() > 1) return 1;
		return data.front() >= abs(a);
	}
	inline bool operator >= (const unsigned &a) {
		if (sign) return 0;
		if (data.size() > 1) return 1;
		return data.front() >= a;
	}
	inline bool operator < (const int &a) {
		if (sign != (a < 0)) return sign;
		if (sign) {
			if (data.size() > 1) return 1;
			return data.front() > abs(a);
		}
		if (data.size() > 1) return 0;
		return data.front() < abs(a);
	}
	inline bool operator < (const unsigned &a) {
		if (sign) return 1;
		if (data.size() > 1) return 0;
		return data.front() < a;
	}
	inline bool operator <= (const int &a) {
		if (sign != (a < 0)) return sign;
		if (sign) {
			if (data.size() > 1) return 1;
			return data.front() >= abs(a);
		}
		if (data.size() > 1) return 0;
		return data.front() <= abs(a);
	}
	inline bool operator <= (const unsigned &a) {
		if (sign) return 1;
		if (data.size() > 1) return 0;
		return data.front() <= a;
	}

	inline integer& operator ++ () {
		int i = 0;
		if (sign) {
			if (data.size() == 1 && data.front() == 1) {
				data.clear(); data.push_back(0); sign = 0; return *this;
			}
			data.front()--;
			while (data.at(i) == -1) {
				i++;
				data.at(i)--;
			}
			return *this;
		}
		data.front()++;
		while (!data.at(i)) {
			if (i >= data.size() - 1) {
				data.push_back(1);
				return *this;
			}
			i++, data.at(i)++;
		}
		return *this;
	}
	inline integer operator ++ (int) {
		int i = 0;
		if (sign) {
			if (data.size() == 1 && data.front() == 1) {
				data.clear(); data.push_back(0); sign = 0; return *this;
			}
			data.front()--;
			while (data.at(i) == -1) {
				i++;
				data.at(i)--;
			}
			return *this;
		}
		data.front()++;
		while (!data.at(i)) {
			if (i >= data.size() - 1) {
				data.push_back(1);
				return *this;
			}
			i++, data.at(i)++;
		}
		return *this;
	}
	inline integer& operator -- () {
		int i = 0;
		if (sign) {
			data.front()++;
			while (!data.at(i)) {
				if (i >= data.size() - 1) {
					data.push_back(1);
					return *this;
				}
				i++, data.at(i)++;
			}
			return *this;
		}
		if (data.size() == 1 && data.front() == 1) {
			data.clear(); data.push_back(0); sign = 0; return *this;
		}
		if (data.size() == 1 && data.front() == 0) {
			data.clear(); data.push_back(1); sign = 1; return *this;
		}
		data.front()--;
		while (data.at(i) == -1) {
			i++;
			data.at(i)--;
		}
		return *this;
	}
	inline integer operator -- (int) {
		int i = 0;
		if (sign) {
			data.front()++;
			while (!data.at(i)) {
				if (i >= data.size() - 1) {
					data.push_back(1);
					return *this;
				}
				i++, data.at(i)++;
			}
			return *this;
		}
		if (data.size() == 1 && data.front() == 1) {
			data.clear(); data.push_back(0); sign = 0; return *this;
		}
		if (data.size() == 1 && data.front() == 0) {
			data.clear(); data.push_back(1); sign = 1; return *this;
		}
		data.front()--;
		while (data.at(i) == -1) {
			i++;
			data.at(i)--;
		}
		return *this;
	}

	friend integer operator + (integer a, integer b) {
		unsigned long long c;
		int size = a.data.size() > b.data.size() ? a.data.size() : b.data.size(); size++;
		a.data.resize(size, 0), b.data.resize(size, 0);
		if (a.sign == b.sign) {
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (c) a.data.push_back(c);
			if (!a.data.back()) a.data.pop_back();
			return a;
		}
		else if (a.sign) {
			for (int i = 0; i < size; i++) {
				a.data.at(i) = ~a.data.at(i);
			}
			a.data.at(0)++;
			int i = 0;
			while (!a.data.at(i)) {
				if (i + 1 >= size) break;
				i++, a.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (a.data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					a.data.at(i) = ~a.data.at(i);
				}
				a.data.at(0)++;
				i = 0;
				while (!a.data.at(i)) {
					if (i >= size) break;
					i++, a.data.at(i)++;
				}
			}
			else a.sign = 0;
			while (a.data.back() == 0) {
				a.data.pop_back();
				if (a.data.empty()) {
					a.data.push_back(0);
					a.sign = 0;
					break;
				}
			}
			return a;
		}
		else if (b.sign) {
			for (int i = 0; i < size; i++) {
				b.data.at(i) = ~b.data.at(i);
			}
			b.data.at(0)++;
			int i = 0;
			while (!b.data.at(i)) {
				if (i + 1 >= size) break;
				i++, b.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (a.data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					a.data.at(i) = ~a.data.at(i);
				}
				a.data.at(0)++;
				i = 0;
				while (!a.data.at(i)) {
					if (i >= size) break;
					i++, a.data.at(i)++;
				}
				a.sign = 1;
			}
			while (a.data.back() == 0) {
				a.data.pop_back();
				if (a.data.empty()) {
					a.data.push_back(0);
					a.sign = 0;
					break;
				}
			}
			return a;
		}
	}
	void operator += (integer a) {
		unsigned long long c;
		int size = data.size() > a.data.size() ? data.size() : a.data.size(); size++;
		data.resize(size, 0), a.data.resize(size, 0);
		if (sign == a.sign) {
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (c) data.push_back(c);
			if (!data.back()) data.pop_back();
			return;
		}
		else if (sign) {
			for (int i = 0; i < size; i++) {
				data.at(i) = ~data.at(i);
			}
			data.at(0)++;
			int i = 0;
			while (!data.at(i)) {
				if (i + 1 >= size) break;
				i++, data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					data.at(i) = ~data.at(i);
				}
				data.at(0)++;
				i = 0;
				while (!data.at(i)) {
					if (i >= size) break;
					i++, data.at(i)++;
				}
			}
			else sign = 0;
			while (data.back() == 0) {
				data.pop_back();
				if (data.empty()) {
					data.push_back(0);
					sign = 0;
					break;
				}
			}
			return;
		}
		else if (a.sign) {
			for (int i = 0; i < size; i++) {
				a.data.at(i) = ~a.data.at(i);
			}
			a.data.at(0)++;
			int i = 0;
			while (!a.data.at(i)) {
				if (i + 1 >= size) break;
				i++, a.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					data.at(i) = ~data.at(i);
				}
				data.at(0)++;
				i = 0;
				while (!data.at(i)) {
					if (i >= size) break;
					i++, data.at(i)++;
				}
				sign = 1;
			}
			while (data.back() == 0) {
				data.pop_back();
				if (data.empty()) {
					data.push_back(0);
					sign = 0;
					break;
				}
			}
			return;
		}
	}
	inline friend integer operator - (integer a) {
		if (a.empty()) return a;
		a.sign ^= 1;
		return a;
	}
	friend integer operator - (integer a, integer b) {
		unsigned long long c;
		int size = a.data.size() > b.data.size() ? a.data.size() : b.data.size(); size++;
		a.data.resize(size, 0), b.data.resize(size, 0);
		if (a.sign != b.sign) {
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (c) a.data.push_back(c);
			return a;
		}
		else if (a.sign) {
			for (int i = 0; i < size; i++) {
				a.data.at(i) = ~a.data.at(i);
			}
			a.data.at(0)++;
			int i = 0;
			while (!a.data.at(i)) {
				if (i + 1 >= size) break;
				i++, a.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (a.data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					a.data.at(i) = ~a.data.at(i);
				}
				a.data.at(0)++;
				i = 0;
				while (!a.data.at(i)) {
					if (i >= size) break;
					i++, a.data.at(i)++;
				}
			}
			else a.sign = 0;
			while (a.data.back() == 0) {
				a.data.pop_back();
				if (a.data.empty()) {
					a.data.push_back(0);
					a.sign = 0;
					break;
				}
			}
			return a;
		}
		else {
			for (int i = 0; i < size; i++) {
				b.data.at(i) = ~b.data.at(i);
			}
			b.data.at(0)++;
			int i = 0;
			while (!b.data.at(i)) {
				if (i + 1 >= size) break;
				i++, b.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += a.data.at(i);
				c += b.data.at(i);
				a.data.at(i) = c, c >>= 32;
			}
			if (a.data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					a.data.at(i) = ~a.data.at(i);
				}
				a.data.at(0)++;
				i = 0;
				while (!a.data.at(i)) {
					if (i >= size) break;
					i++, a.data.at(i)++;
				}
				a.sign = 1;
			}
			while (a.data.back() == 0) {
				a.data.pop_back();
				if (a.data.empty()) {
					a.data.push_back(0);
					a.sign = 0;
					break;
				}
			}
			return a;

		}
	}
	void operator -= (integer a) {
		unsigned long long c;
		int size = data.size() > a.data.size() ? data.size() : a.data.size(); size++;
		data.resize(size, 0), a.data.resize(size, 0);
		if (sign != a.sign) {
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (c) data.push_back(c);
			return;
		}
		else if (sign) {
			for (int i = 0; i < size; i++) {
				data.at(i) = ~data.at(i);
			}
			data.at(0)++;
			int i = 0;
			while (!data.at(i)) {
				if (i + 1 >= size) break;
				i++, data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					data.at(i) = ~data.at(i);
				}
				data.at(0)++;
				i = 0;
				while (!data.at(i)) {
					if (i >= size) break;
					i++, data.at(i)++;
				}
			}
			else sign = 0;
			while (data.back() == 0) {
				data.pop_back();
				if (data.empty()) {
					data.push_back(0);
					sign = 0;
					break;
				}
			}
			return;
		}
		else {
			for (int i = 0; i < size; i++) {
				a.data.at(i) = ~a.data.at(i);
			}
			a.data.at(0)++;
			int i = 0;
			while (!a.data.at(i)) {
				if (i + 1 >= size) break;
				i++, a.data.at(i)++;
			}
			c = 0;
			for (int i = 0; i < size; i++) {
				c += data.at(i);
				c += a.data.at(i);
				data.at(i) = c, c >>= 32;
			}
			if (data.back() >> 31) {
				for (int i = 0; i < size; i++) {
					data.at(i) = ~data.at(i);
				}
				data.at(0)++;
				i = 0;
				while (!data.at(i)) {
					if (i >= size) break;
					i++, data.at(i)++;
				}
				sign = 1;
			}
			while (data.back() == 0) {
				data.pop_back();
				if (data.empty()) {
					data.push_back(0);
					sign = 0;
					break;
				}
			}
			return;
		}
	}
	friend integer operator * (integer a, integer b) {
		// Better Algorithm: https://www.cnblogs.com/zuoxiaolong/p/computer10.html
		unsigned long long g;
		integer c;
		if (a.empty() || b.empty()) return c;
		c.data.clear(); c.sign = a.sign ^ b.sign;
		if (a.data.size() == 1 && b.data.size() == 1) {
			g = a.data.front(); g *= b.data.front();
			c.data.push_back(g); g >>= 32;
			if (g) c.data.push_back(g);
			return c;
		}
		c.data.resize(a.data.size() + b.data.size());
		int n;
		for (int i = 0; i < b.data.size(); i++) {
			g = 0;
			for (int j = 0; j < a.data.size(); j++) {
				g += (unsigned long long)(a.data.at(j))*(unsigned long long)(b.data.at(i));
				g += c.data.at(i + j);
				c.data.at(i + j) = g, g >>= 32;
				n = i + j + 1;
				while (g) {
					g += c.data.at(n), c.data.at(n) = g, g >>= 32, n++;
				}
			}
		}
		while (!c.data.back()) c.data.pop_back();
		return c;
	}
	inline void operator *= (const integer &a) {
		*this = (*this)*a;
	}
	friend integer operator / (integer a, integer b) {
		integer c;
		c.sign = a.sign ^ a.sign, a.sign = b.sign = 0;
		if (a < b) {
			c.sign = 0; return c;
		}
		//if (b.empty());
		unsigned long long m, n;
		if (a.data.size() <= 2 && b.data.size() <= 2) {
			a.data.resize(2, 0), b.data.resize(2, 0);
			m = a.data.back(), m <<= 32, m |= a.data.front();
			n = b.data.back(), n <<= 32, n |= b.data.front();
			m /= n;
			c.data.front() = m;
			m >>= 32;
			if (m) c.data.push_back(m);
			return c;
		}
		int t = a.data.size() - b.data.size();
		c.data.resize(t + 1);
		unsigned long long D = b.data.back();
		if (b.data.size() == 1) {
			m = 0;
			for (int i = t; i >= 0; i--) {
				m <<= 32; m |= a.data.at(i);
				c.data.at(i) = m / D;
				m %= D;
			}
			if (!c.data.back()) c.data.pop_back();
			return c;
		}
		D <<= 32; D |= b.data.at(b.data.size() - 2);
		if (a.data.size() == b.data.size()) {
			m = a.data.back(), m <<= 32, m |= a.data.at(a.data.size() - 2);
			n = b.data.back(), n <<= 32, n |= b.data.at(b.data.size() - 2);
			c.data.front() = m / n;
			if (n >= 0x8000000000000000) return c;
			if (c.data.front() != (m / (n + 1))) {

			}
			return c;
		}
		for (int i = t; i >= 0; i--) {
			// floor[(a1*k+a2)/(b1*k+b2+1)] ~ floor[(a1*k+a2)/(b1*k+b2)]    k=2^32, limit of unsigned int

		}

	}
	inline void operator /= (integer a) {
		*this = *this / a;
	}
	//friend integer operator % (integer a, integer b);
	//void operator %= (integer a);
	//integer pow(integer a, integer b);
	//integer pow(integer a, integer b, integer c);

	//inline friend integer operator + (integer a, const int &b);
	inline friend integer operator + (integer a, const unsigned &b) {
		unsigned long long c = b;
		int i = 0;
		if (a.sign) {
			while (c) {

			}
		}
		while (c) {
			c += a.data.at(i);
			a.data.at(i) = c, c >>= 32;
			if (c && i + 1 >= a.data.size()) {
				a.data.push_back(c);
				break;
			}
			if (!c) break;
			i++;
		}
		return a;
	}
	//inline void operator += (const int &a);
	//inline void operator += (const unsigned &a);
	//inline friend integer operator - (integer a, const int &b);
	//inline friend integer operator - (integer a, const int &b);
	//inline friend integer operator - (integer a, const unsigned &b);
	//inline void operator -= (const int &b);
	//inline void operator -= (const unsigned &b);
	//inline friend integer operator * (integer a, const int &b);
	//inline friend integer operator * (integer a, const unsigned &b);
	//inline void operator *= (const int &b);
	//inline void operator *= (const unsigned &b);
	//inline friend integer operator / (integer a, const int &b);
	//inline friend integer operator / (integer a, const unsigned &b);
	//inline void operator /= (const int &b);
	//inline void operator /= (const unsigned &b);
	//inline friend integer operator % (integer a, const int &b);
	//inline friend integer operator % (integer a, const unsigned &b);
	//inline void operator %= (const int &b);
	//inline void operator %= (const unsigned &b);
	explicit operator unsigned() const {
		return data.front();
	}
	explicit operator int() const {
		if (sign) return -int(data.front());
		return data.front();
	}
	friend inline integer pow(integer a, unsigned b) {
		integer c = 1;
		while (b) {
			if (b & 1) c *= a;
			a *= a;
			b >>= 1;
		}
		return c;
	}
	//inline integer pow(integer a, unsigned b, integer c);
	//inline integer pow(integer a, unsigned b, int c);

	//inline friend integer operator & (integer a, integer b);
	//inline integer operator &= (integer a);
	//inline friend integer operator | (integer a, integer b);
	//inline integer operator |= (integer a);
	//inline friend integer operator ^ (integer a, integer b);
	//inline integer operator ^= (integer a);
	//inline friend integer operator ~ (integer a);
	//inline friend integer operator << (integer a, const unsigned &b);
	//inline integer operator <<= (const unsigned &a);
	//inline friend integer operator >> (integer a, const unsigned &b);
	//inline integer operator >>= (const unsigned &a);
	//inline friend integer operator & (integer a, const unsigned &b);
	//inline integer operator &= (const unsigned &a);
	//inline friend integer operator | (integer a, const unsigned &b);
	//inline integer operator |= (const unsigned &a);
	//inline friend integer operator ^ (integer a, const unsigned &b);
	//inline integer operator ^= (const unsigned &a);

	//bool prime();
	//friend integer gcf(integer a, integer b);
	//friend integer lcm(integer a, integer b);
	//friend integer gcd(integer a, integer b, integer c);
	//inline bool divisible(const integer &divisor);
	//inline bool congruence(const integer &divisor, const integer remainder);

	friend class fraction;

	/* For Debug */
	friend inline integer RandomInteger_NormalDistribution();
	friend inline integer RandomInteger_NearNode();
	friend int main();
};

#include "D:\Coding\AboutMath\Calculator\Calculator\Define.h"
double erfinv(double x) {
	double t, n;
	n = log(1 - x * x);
	t = 2 / (_const_pi*0.147) + 0.5 * n;
	if (x > 0) return sqrt(-t + sqrt(t*t - (1 / (0.147) * n)));
	return -sqrt(-t + sqrt(t*t - (1 / (0.147) * n)));
}
double randnor(double median, double variance) {
	long long k = 0;
	for (int i = 0; i < 64; i++) {
		k |= rand() & 1;
		k <<= 1;
	}
	double d = double(k) / 9223372036854775810.0;	// d in (-1,1)
	return erfinv(d)*sqrt(2)*variance + median;
}
#include <chrono>
class _MyT_RANDOM {
public:
	_MyT_RANDOM() {
		chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds>(
			chrono::system_clock::now().time_since_epoch()
			);
		srand((unsigned)this ^ (unsigned)(ms.count()));
	}
};
_MyT_RANDOM _INIT_RANDOM;
inline unsigned Rand() {
	unsigned k = 0;
	for (int i = 0; i < 32; i++) {
		k <<= 1; k |= rand() & 1;
	}
	return k;
}
#define RandomInteger() RandomInteger_NormalDistribution()
inline integer RandomInteger_NormalDistribution() {
	integer a;
	int s = ceil(64 * erfinv((Rand() - 1) / 4294967296.0) + 1);
	a.data.resize(s / 32 + 1, 0);
	for (int i = 0; i < a.data.size() - 1; i++) {
		for (int j = 0; j < 32; j++) {
			a.data.at(i) |= rand() & 1 << j;
		}
	}
	s %= 32;
	for (int i = 0; i < s; i++) {
		a.data.back() |= rand() & 1 << i;
	}
	a.sign = rand() & 1;
	return a;
	/*
	+-----------+-------------+
	| data.size | probability |
	+-----------+-------------+
	|     1     |  52.05%     |
	|     2     |  32.22%     |
	|     3     |  12.34%     |
	|     4     |  2.922%     |
	|     5     |  0.427%     |
	|     6     |  0.0385%    |
	|     7     |  0.00214%   |
	|     8     |  0.000073%  |
	+-----------+-------------+
	*/
}
inline integer RandomInteger_NearNode() {
	integer a;
	int s = ceil(2 * erfinv((Rand() + 1.0) / 4294967300.0));
	a.data.resize(s, 0);
	for (int i = 0; i < s; i++) {
		a.data.at(i) = unsigned(floor(4 * erfinv((Rand() - 1) / 4294967296.0) + 1));
		if (rand() & 1) {
			a.data.at(i) = ~a.data.at(i); a.data.at(i)++;
		}
	}
	a.sign = rand() & 1;
	return a;
}

#endif



#ifndef _INC_FRACTION

#define _INC_FRACTION

#ifdef LargeInteger
class fraction {
	integer den, num;
public:
	fraction() :den(1), num(0) {};
	template<typename T> fraction(const T& numerator) :num(numerator) {}
	template<typename T1, typename T2> fraction(const T1 &numerator, const T2 &denominator) {
		den = denominator, num = numerator; simplify();
	}
	fraction(const fraction &other) :num(other.num), den(other.den) { simplify(); }
	fraction(const string &a) {
		*this = a;
	}
	template<typename T> inline fraction& operator = (const T& number) {
		den = 1, num = number;
		return *this;
	}
	inline fraction& operator = (const fraction &other) {
		num = other.num, den = other.den; simplify();
		return *this;
	}
	fraction& operator = (string a) {
		bool sign = 0;
		if (a[0] == '-') sign = 1, a.erase(0, 1);
		if (a[a.find('/', 0) + 1] == '-' && a.find('/', 0) != -1) sign ^= 1;
		if (a.find('/', 0) == -1) {
			if (a.find('.', 0) != -1) {
				for (int i = 0; i < a.size(); i++) {
					if ((a[i] < 47 || a[i] > 57) && a[i] != '.') a.erase(i, 1);
				}
				for (int i = a.find('.', 0) + 1; i < a.size(); i++) {
					if (a[i] == '.') a.erase(i, 1), i--;
				}
				den = 1;
				for (int i = a.size() - a.find('.', 0) - 1; i > 0; i--) {
					den *= 10;
				}
				a.erase(a.find('.', 0), 1);
				num = a;
				if (sign) num = -num;
				simplify();
				return *this;
			}
			for (int i = 0; i < a.size(); i++) {
				if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
			}
			num = a, den = 1;
			if (sign) num = -num;
			simplify();
			return *this;
		}
		for (int i = 0; i < a.size(); i++) {
			if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
		}
		num = a.substr(0, a.find('/', 0));
		den = a.substr(a.find('/', 0) + 1, a.size() - a.find('/', 0));
		if (sign) num.minus();
		simplify();
		return *this;
	}
	inline friend istream& operator >> (istream& is, fraction &a)
	{
		string n; is >> n;
		a = n;
		return is;
	}
	inline friend ostream& operator << (ostream& os, fraction a)
	{
		if (a.den.empty()) {
			if (a.num > 0) os << "#INF";
			else if (a.num < 0) os << "#-INF";
			else if (a.num.empty()) os << "#NAF";
			return os;
			// "#INF" will become "#NAN" after any calculation. 
		}
		a.simplify();
		if (a.num.empty()) {
			os << "0";
			return os;
		}
		if (a.den == 1) {
			os << a.num;
			return os;
		}
		os << a.num << "/" << a.den;
		return os;
	}

	inline void simplify() {
		if (den.empty()) {
			num.clear();
			return;
		}
		if (num.empty()) {
			den = 1; return;
		}
		//integer m = gcf(num, den);
		//if (m == 1) return;
		//den /= m; num /= m;
	}

	inline bool proper() {
		return abs(num) < abs(den);
	}
	inline bool imporper() {
		return abs(num) >= abs(den);
	}
	inline bool whole() {
		simplify(); return den == 1;
	}
	inline integer denominator() {
		return den;
	}
	inline integer numerator() {
		return num;
	}

};
#else

unsigned _Myt_Fraction_GCF(unsigned a, unsigned b)
{
	if (a == 0 || b == 0) return 0;
	unsigned c = -1;
	while (c != 0) c = a % b, a = b, b = c;
	return a;
}
unsigned _Myt_Fraction_LCM(unsigned a, unsigned b)
{
	if (a == 0 || b == 0) return 0;
	return b / _Myt_Fraction_GCF(a, b) * a;
}
unsigned _Myt_Fraction_String_to_Num(string a) {
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
				for (int i = 0; i < a.size(); i++) {
					if ((a[i] < 47 || a[i] > 57) && a[i] != '.') a.erase(i, 1);
				}
				for (int i = a.find('.', 0) + 1; i < a.size(); i++) {
					if (a[i] == '.') a.erase(i, 1), i--;
				}
				den = 1;
				for (int i = a.size() - a.find('.', 0) - 1; i > 0; i--) {
					den *= 10;
				}
				a.erase(a.find('.', 0), 1);
				num = _Myt_Fraction_String_to_Num(a);
				if (sign) num = -num;
				simplify();
				return;
			}
			for (int i = 0; i < a.size(); i++) {
				if (a[i] < 47 || a[i] > 57) a.erase(i, 1);
			}
			num = _Myt_Fraction_String_to_Num(a), den = 1;
			if (sign) num = -num;
			simplify();
			return;
		}
		for (int i = 0; i < a.size(); i++) {
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
		return b * a.den < a.num;
	}
	inline friend bool operator < (fraction a, int b) {
		return b * a.den > a.num;
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

#endif


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
	matrix(matrix<T> &M) {
		h = M.h, w = M.w;
		if (h == 0 || w == 0) { data = 0; return; }
		data = new T[h * w];
		for (int i = h * w; i > 0; i--) {
			*data = *M.data;
			data++, M.data++;
		}
		data -= h * w, M.data -= h * w;
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
	matrix<T>& operator = (matrix<T> &M) {
		if (this->h != M.h || this->w != M.w) {
			this->h = M.h, this->w = M.w;
			if (data != 0) delete data;
			data = new T[h*w];
		}
		for (int i = h * w; i > 0; i--) {
			*data = *M.data;
			data++, M.data++;
		}
		data -= h * w, M.data -= h * w;
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
	inline const unsigned height() { return this->h; }
	inline const unsigned width() { return this->w; }
	inline const unsigned size() { return h * w; }
	inline const unsigned length() { return h * w; }
	friend ostream& operator << (ostream& os, const matrix<T> &M)
	{
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
		return os;
	}

	inline bool empty() {
		return h == 0 || w == 0 || data == 0;
	}
	inline bool square() {
		return h == w && h != 0;
	}
	inline bool row() {
		return h == 1 && w != 0;
	}
	inline bool column() {
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
			for (int j = 0; j < i; j++) {
				if ((*this)[i][j] != (*this)[j][i]) return 0;
			}
		}
	}
	bool skew_symmetric() {
		if (h != w) return 0;
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < i; j++) {
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
		if (h != w) return 0;
		if (h == 0) return data == 0 ? 0 : 1;
		if (h == 1) return *data;
		if (h == 2) return (*data)*(*(data + 3)) - (*(data + 1)*(*(data + 2)));
		if (h == 3) return (*data)*((*(data + 4))*(*(data + 8)) - (*(data + 5)*(*(data + 7))))
			- (*(data + 1))*((*(data + 3))*(*(data + 8)) - (*(data + 5)*(*(data + 6))))
			+ (*(data + 2))*((*(data + 3))*(*(data + 7)) - (*(data + 4)*(*(data + 6))));
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
						//cout << "MultiAdd E: " << endl << *this << endl;
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
					//cout << "MultiAdd M: " << endl << *this << endl;
				}
			}
		}
		//cout << *this << endl;
		p = A.data; c = 1;
		for (int i = 0; i < h; i++) {
			if (*p == 0) return 0;
			c *= *p, p += w + 1;
		}
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

	virtual vector<T> solve(vector<T> r) {
		return vector<T>();
	}


#ifdef OLD_SOURCE
	void elimination_old() {
		if (h == 0 || w == 0 || h == 1) return;
		if (w == 1) {
			for (int i = 0; i < h; i++) {
				*data = 0, data++;
			}
			data -= h;
			*data = 1;
			return;
		}
		T *p, *q, c;
		for (int i = 0, j = 0; i < w; i++) {
			while ((*this)[j][i] == 0) {	// begin with 0
				for (int k = j + 1; k < h; k++) {
					if ((*this)[k][i] != 0) {
						// Interchange(j, k)
						p = data + j * w + i, q = data + k * w + i;
						for (int n = i; n < w; n++) {
							c = *p, *p = *q, *q = c;
							p++, q++;
						}
						p = q = 0;
						//cout << "Interchange the " << (j + 1) << "th and the " << (k + 1) << "th row. \n" << (*this) << endl;
						break;
					}
					if (k + 1 >= h) {	// whole colume is 0
						i++; break;
					}
				}
				if (i + 1 >= w) break;
				if (j + 1 >= h) break;
			}
			for (int k = j + 1; k < h; k++) {
				if ((*this)[k][i] != 0) {
					// Multiple j and Add to k
					c = -(*this)[k][i] / (*this)[j][i];
					p = data + j * w + i, q = data + k * w + i;
					for (int n = i; n < w; n++) {
						*q += c * (*p);
						q++, p++;
					}
					p = q = 0;
					//cout << "Add " << c << " times the " << (j + 1) << "th row to the " << (k + 1) << "th row. \n" << (*this) << endl;
				}
			}
			j++;
			if (j >= h) break;
		}
		for (int i = h - 1; i > 0; i--) {
			p = data + i * w;
			int k;
			for (k = 0; k < w; k++) {
				if (*p != 0) break;
				p++;
			}
			if (k != w) {
				for (int j = i - 1; j >= 0; j--) {
					if ((*this)[j][k] != 0) {
						c = -((*this)[j][k] / (*this)[i][k]);
						p = data + i * w + k, q = data + j * w + k;
						for (int n = k; n < w; n++) {
							*q += c * (*p);
							q++, p++;
						}
						p = q = 0;
						//cout << "Add " << c << " times the " << (i + 1) << "th row to the " << (j + 1) << "th row. \n" << (*this) << endl;
					}
				}
			}
		}
		p = data;
		for (int i = 0; i < h; i++) {
			c = 0;
			for (int j = 0; j < w; j++) {
				if (*p != 0) {
					if (c == 0) c = 1 / (*p), *p = 1;
					else (*p) *= c;

				}
				p++;
			}
		}
		p = 0;
		//cout << "Final Matrix: \n" << (*this) << endl << endl;
	}

	/*
	Main()
	cout << "Matrix A: \n" << use << endl;
	cout << "Exponent calculate with Taylor Series (160 terms): \n" << exp_s(use) << endl;
	cout << "Sine calculate with Taylor Series (160 terms): \n" << sin_s(use) << endl;
	cout << "Cosine calculate with Taylor Series (160 terms): \n" << cos_s(use) << endl;
	cout << "sin^2(A)+cos^2(A): \n" << (pow(sin_s(use), 2) + pow(cos_s(use), 2)) << endl;
	cout << Exp(use) << endl;
	cout << Sin(use) << endl;
*/


#endif
#ifdef REC_SOURCE
	friend matrix<T> Exp(matrix<T> A) {
		matrix<T> R;
		if (A.h != A.w) return R;
		if (A.h == 2) {
			R.h = R.w = 2;
			R.data = new T[4];
			T a = *A.data, b = *(A.data + 1), c = *(A.data + 2), d = *(A.data + 3);
			T D = (a - d)*(a - d) + 4 * b*c;
			T e = exp((a + d) / 2);
			if (D > 0) {
				D = sqrt(D);
				*R.data = (D*cosh(D / 2) + (a - d)*sinh(D / 2))*e / D;
				*(R.data + 1) = 2 * b*e*sinh(D / 2) / D;
				*(R.data + 2) = 2 * c*e*sinh(D / 2) / D;
				*(R.data + 3) = (D*cosh(D / 2) + (d - a)*sinh(D / 2))*e / D;
				return R;
			}
			else if (D == 0) {
				*R.data = (1 + (a - d) / 2)*e;
				*(R.data + 1) = b * e;
				*(R.data + 2) = c * e;
				*(R.data + 3) = (1 - (a - d) / 2)*e;
				return R;
			}
			else {
				// sinh(bi)=sin(b)i, cosh(bi)=cos(b)
				D = sqrt(-D);
				*R.data = (D*cos(D / 2) + (a - d)*sin(D / 2))*e / D;
				*(R.data + 1) = 2 * b*e*sin(D / 2) / D;
				*(R.data + 2) = 2 * c*e*sin(D / 2) / D;
				*(R.data + 3) = (D*cos(D / 2) + (d - a)*sin(D / 2))*e / D;
				return R;
			}
		}
		return R;
	}
	friend matrix<T> Sin(matrix<T> A) {
		matrix<T> R;
		if (A.h != A.w) return R;
		if (A.h == 2) {
			R.h = R.w = 2;
			R.data = new T[4];
			T a = *A.data, b = *(A.data + 1), c = *(A.data + 2), d = *(A.data + 3);		// I
			T D = (a - d)*(d - a) - 4 * b*c;	// R
			T er = cos((a + d) / 2), ei = sin((a + d) / 2);		// C
			if (D > 0) {
				D = sqrt(D);
				*R.data = (D*cosh(D / 2) + (a - d)*sinh(D / 2))*e / D;
				*(R.data + 1) = 2 * b*e*sinh(D / 2) / D;
				*(R.data + 2) = 2 * c*e*sinh(D / 2) / D;
				*(R.data + 3) = (D*cosh(D / 2) + (d - a)*sinh(D / 2))*e / D;
				return R;
			}
			else if (D == 0) {
				*R.data = (1 + (a - d) / 2)*e;
				*(R.data + 1) = b * e;
				*(R.data + 2) = c * e;
				*(R.data + 3) = (1 - (a - d) / 2)*e;
				return R;
			}
			else {
				// sinh(bi)=sin(b)i, cosh(bi)=cos(b)
				D = sqrt(-D);
				*R.data = (D*cos(D / 2) + (a - d)*sin(D / 2))*e / D;
				*(R.data + 1) = 2 * b*e*sin(D / 2) / D;
				*(R.data + 2) = 2 * c*e*sin(D / 2) / D;
				*(R.data + 3) = (D*cos(D / 2) + (d - a)*sin(D / 2))*e / D;
				return R;
			}
		}
		return R;
	}
	friend matrix<T> exp_s(matrix<T> A) {
		if (!A.sqr()) return matrix();
		matrix<T> M(A.h, Matrix_type::_M_Identity);
		matrix<T> F = M;
		matrix<T> R(A.h, A.w);
		double fact = 1;
		for (int i = 0; i < 160; i++) {
			R += F;
			//cout << i << " " << fact << " " << determinant(M) << endl << R << endl;
			fact *= i + 1;
			M *= A;
			F = M / fact;
			if (F.zero()) break;
		}
		return R;
	}
	friend matrix<T> sin_s(matrix<T> A) {
		if (!A.sqr()) return matrix();
		matrix<T> M(A.h, Matrix_type::_M_Identity);
		matrix<T> F = M;
		matrix<T> R(A.h, A.w);
		double fact = 1;
		for (int i = 0; i < 160; i++) {
			if (i & 1) {
				if ((i / 2) & 1) R -= F;
				else R += F;
			}
			//cout << i << " " << fact << " " << determinant(M) << endl << R << endl;
			fact *= i + 1;
			M *= A;
			F = M / fact;
			if (M.zero()) break;
		}
		return R;
	}
	friend matrix<T> cos_s(matrix<T> A) {
		if (!A.sqr()) return matrix();
		matrix<T> M(A.h, Matrix_type::_M_Identity);
		matrix<T> F = M;
		matrix<T> R(A.h, A.w);
		double fact = 1;
		for (int i = 0; i < 160; i++) {
			if (!(i & 1)) {
				if ((i / 2) & 1) R -= F;
				else R += F;
			}
			//cout << i << " " << fact << " " << determinant(M) << endl << R << endl;
			fact *= i + 1;
			M *= A;
			F = M / fact;
			if (F.zero()) break;
		}
		return R;
	}
#endif

};

template<typename T> class Vandermonde :public matrix<T> {
	vector<T> v;
public:
	inline Vandermonde() :h(0), w(0), data(0) {}
	Vandermonde(Vandermonde<T> &M) {
		h = M.h, w = M.w;
		if (h == 0 || w == 0) { data = 0; return; }
		data = new T[h * w];
		for (int i = h * w; i > 0; i--) {
			*data = *M.data;
			data++, M.data++;
		}
		data -= h * w, M.data -= h * w;
		if (!v.empty()) v.clear();
		for (int i = 0; i < h; i++) {
			v.push_back(M.v.at(i));
		}
	}
	Vandermonde(const vector<T> &q) {
		int ps = q.size();
		h = w = ps;
		data = new T[h*w];
		T n;
		for (int i = 0; i < h; i++) {
			n = q.at(i);
			for (int j = 0; j < w; j++) {
				*data = pow(n, j);
				data++;
			}
		}
		data -= h * w;
		v = q;
	}
	Vandermonde(const vector<T> &q, Matrix_type::_M_type type) {
		int ps = q.size();
		h = w = ps;
		data = new T[h*w];
		T n;
		v = q;
		if (type.type & 1 || type.type & 2) {
			for (int i = 0; i < h; i++) {
				n = q.at(i);
				for (int j = 0; j < w; j++) {
					*data = pow(n, j);
					data++;
				}
			}
			data -= h * w;
			return;
		}
		else if (type.type & 4) {
			for (int i = 0; i < h; i++) {
				for (int j = 0; j < w; j++) {
					*data = pow(q.at(j), i);
					data++;
				}
			}
			data -= h * w;
			return;
		}
		else {
			for (int i = h * w; i > 0; i--) {
				*data = 0; data++;
			}
			data -= h * w;
			v.clear();
			return;
		}
	}
	Vandermonde(const initializer_list<T> &q) {
		int ps = q.size();
		h = w = ps;
		int ds = h * w;
		data = new T[ds];
		v.clear();
		const T* p = q.begin();
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				*data = pow(*p, j);
				data++;
			}
			v.push_back(*p);
			p++;
		}
		data -= ds;
	}
	Vandermonde(const initializer_list<T> &q, Matrix_type::_M_type type) {
		int ps = q.size();
		h = w = ps;
		int ds = h * w;
		data = new T[ds];
		v.clear();
		const T* p = q.begin();
		if (type.type & 1 || type.type & 2) {
			for (int i = 0; i < h; i++) {
				for (int j = 0; j < w; j++) {
					*data = pow(*p, j);
					data++;
				}
				v.push_back(*p);
				p++;
			}
			data -= ds;
			return;
		}
		else if (type.type & 4) {
			for (int i = 0; i < h; i++) {
				for (int j = 0; j < w; j++) {
					*data = pow(*p, i);
					data++, p++;
				}
				p -= w;
			}
			do {
				v.push_back(*p);
			} while (++p < q.end());
			data -= ds;
			return;
		}
		else {
			for (int i = ds; i > 0; i--) {
				*data = 0; data++;
			}
			data -= ds;
			return;
		}
	}
	~Vandermonde() {
		for (int i = 0; i < v.size(); i++) {
			//cout << v.at(i) << " ";
		}
		h = w = 0;
		delete data; data = 0;
		v.clear();
	}
	inline friend ostream& operator << (ostream& os, const Vandermonde<T> &M) {
		os << *((matrix<T>*)&M);
		return os;
	}

	template<typename t> inline friend matrix<T> operator * (const t &b, Vandermonde<T> a) { return b * (*((matrix<T>*)&a)); }
	inline friend matrix<T> pow(Vandermonde<T> a, int b) { return pow(*((matrix<T>*)&a), b); }
	inline friend T det(Vandermonde &A) { return A.Det(); }
	inline friend T determinant(Vandermonde &A) { return A.Determinant(); }
	inline friend matrix<T> adj(Vandermonde &A) { return A.adjugate(); }
	inline friend unsigned R(Vandermonde &A) { return A.rank(); }
	inline friend unsigned Rank(Vandermonde A) { return Rank(*((matrix<T>*)&A)); }

	T Det() {
		if (v.empty()) return 0;
		T r = 1;
		for (int i = h - 1; i >= 0; i--) {
			for (int j = i + 1; j < h; j++) {
				r *= v.at(j) - v.at(i);
			}
		}
		return r;
	}
	unsigned rank() {
		bool sw = 0; int sz;
		vector<T> r = v;
		do {
			sw = 0; sz = r.size();
			for (int i = 1; i < sz; i++) {
				if (r.at(i - 1) > r.at(i)) swap(r.at(i - 1), r.at(i)), sw = 1;
				else if (r.at(i - 1) == r.at(i)) r.erase(r.begin() + i), sz--, i--;
			}
		} while (sw);
		return r.size();
	}

	vector<T> solve(vector<T> r) {
		T d = 1, dv;
		for (int i = h - 1; i >= 0; i--) {
			for (int j = i + 1; j < h; j++) {
				d *= v.at(j) - v.at(i);
			}
		}
		vector<T> R;
		if (d == 0) return R;
		for (int i = 1; i <= w; i++) {
			dv = 0;
			for (int j = 1; j <= h; j++) {
				if ((i + j) & 1) dv -= r.at(j - 1) * this->submatrix(j, i).Determinant();
				else dv += r.at(j - 1) * this->submatrix(j, i).Determinant();
			}
			R.push_back(dv / d);
		}
		return R;
	}

};

template<typename T> class Vector {
	unsigned d;
	T* data;
public:
	Vector() :d(0), data(0) {}
	Vector(unsigned demension) :d(demension) {
		data = new T[d];
		for (int i = 0; i < d; i++) data[i] = 0;
	}
	Vector(unsigned demension, initializer_list<T> q) :d(demension) {
		data = new T[d];
		int size = q.size();
		if (d < size()) size = d;
		const T *p = q.begin();
		int i;
		for (i = 0; i < size; i++) data[i] = *p, p++;
		for (; i < d; i++) data[i] = 0;
	}
	Vector(const initializer_list<T> &q) :d(q.size()) {
		data = new T[d];
		const T *p = q.begin();
		for (int i = 0; i < d; i++) data[i] = *p, p++;
	}
	Vector(const Vector<T> &v) {
		d = v.d;
		data = new T[d];
		for (int i = 0; i < d; i++) data[i] = v.data[i];
	}
	Vector(const T &a, const T &b) {
		d = 2; data = new T[2];
		data[0] = a, data[1] = b;
	}
	Vector(const T &a, const T &b, const T &c) {
		d = 3; data = new T[3];
		data[0] = a, data[1] = b, data[2] = c;
	}
	Vector<T>& operator = (const Vector<T> &v) {
		if (data != 0) delete data;
		d = v.d;
		data = new T[d];
		for (int i = 0; i < d; i++) data[i] = v.data[i];
		return *this;
	}
	Vector<T>& operator = (const initializer_list<T> q) {
		if (data != 0) delete data;
		d = q.size();
		data = new T[d];
		const T *p = q.begin();
		for (int i = 0; i < d; i++) data[i] = *p, p++;
		return *this;
	}
	~Vector() {
		if (data != 0) delete data, data = 0;
		d = 0;
	}
	friend ostream& operator << (ostream& os, const Vector<T> &M) {
		if (M.data == 0) { os << "#NAV"; return os; }
		if (M.d == 0) { os << "<empty>"; return os; }
		os << "(";
		for (int i = 0; i < M.d; i++) os << M.data[i] << ",";
		os << "\b)";
		return os;
	}

	inline T& x() { return data[0]; }
	inline T& y() { return data[1]; }
	inline T& z() { return data[2]; }
	inline T& at(int i) { return data[i - 1]; }
	inline constexpr int dimension() { return d; }
	inline T Mod() {
		T r = 0;
		for (int i = 0; i < d; i++) r += data[i] * data[i];
		return sqrt(r);
	}
	friend inline T mod(const Vector<T> &V) {
		T r = 0;
		for (int i = 0; i < V.d; i++) r += V.data[i] * V.data[i];
		return sqrt(r);
	}

	inline Vector<T> operator + (Vector<T> v) {
		if (v.d != d) return Vector();
		for (int i = 0; i < d; i++) v.data[i] += data[i];
		return v;
	}
	inline void operator += (const Vector<T> &v) {
		if (v.d != d) { delete data; data = 0; d = 0; return; }
		for (int i = 0; i < d; i++) data[i] += v.data[i];
	}
	inline Vector<T> operator - () {
		Vector<T> v;
		v.d = d; v.data = new T[d];
		for (int i = 0; i < d; i++) v.data[i] = -data[i];
		return v;
	}
	inline Vector<T> operator - (Vector<T> v) {
		if (v.d != d) return Vector();
		for (int i = 0; i < d; i++) v.data[i] = data[i] - v.data[i];
		return v;
	}
	inline void operator -= (const Vector<T> &v) {
		if (v.d != d) { ~Vector(); return; }
		for (int i = 0; i < d; i++) data[i] -= v.data[i];
	}
	inline Vector<T> operator * (const T &k) {
		Vector<T> v;
		v.d = d; v.data = new T[d];
		for (int i = 0; i < d; i++) v.data[i] = data[i] * k;
		return v;
	}
	inline friend Vector<T> operator * (const T &k, Vector<T> v) {
		for (int i = 0; i < v.d; i++) {
			v.data[i] *= k;
		}
		return v;
	}
	inline Vector<T> operator / (const T &k) {
		Vector<T> v;
		v.d = d; v.data = new T[d];
		for (int i = 0; i < d; i++) v.data[i] = data[i] / k;
		return v;
	}
	inline void operator *= (const T &k) {
		for (int i = 0; i < d; i++) data[i] *= k;
	}
	inline void operator /= (const T &k) {
		for (int i = 0; i < d; i++) data[i] /= k;
	}

	/* Dot Product */
	inline T operator * (const Vector<T> &v) {
		T r = 0;
		for (int i = min(d, v.d) - 1, i >= 0; i--) {
			r += data[i] * v.data[i];
		}
		return r;
	}
	friend T dot(const Vector<T> &u, const Vector<T> &v) {
		T r = 0;
		for (int i = min(d, v.d) - 1, i >= 0; i--) {
			r += u.data[i] * v.data[i];
		}
		return r;
	}

	/* Cross Product */
	inline void operator *= (const Vector<T> &v) {
		Vector<T> r(3);
		if (d >= 3 && v.d >= 3) {
			r.data[0] = data[1] * v.data[2] - data[2] * v.data[1];
			r.data[1] = data[2] * v.data[0] - data[0] * v.data[2];
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else if (d == 2 && v.d == 2) {
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else if (d >= 3 && v.d == 2) {
			r.data[0] = -data[2] * v.data[1];
			r.data[1] = data[2] * v.data[0];
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else if (d == 2 && v.d >= 3) {
			r.data[0] = data[1] * v.data[2];
			r.data[1] = -data[0] * v.data[2];
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else {
			~Vector(); return;
		}
		if (d != 3) {
			delete data; data = new T[3];
		}
		data[0] = r.data[0], data[1] = r.data[1], data[2] = r.data[2];
		return;
	}
	friend T cross(const Vector<T> &u, const Vector<T> &v) {
		if (u.d != v.d) return Vector<T>();
		Vector<T> r(3);
		if (d >= 3 && v.d >= 3) {
			r.data[0] = data[1] * v.data[2] - data[2] * v.data[1];
			r.data[1] = data[2] * v.data[0] - data[0] * v.data[2];
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else if (d == 2 && v.d == 2) {
			r.data[2] = data[0] * v.data[1] - data[1] - v.data[0];
		}
		else return Vector<T>();
		return r;
	}

	/* Mixed Product */
	/*inline T operator [] (const Vector<T> &u, const Vector<T> &v) {
		T r = 0;
		if (u.d != 2 || v.d != 2) return r;
		r = (*u.data)*v.data[1] - u.data[1] * (*v.data);
		return r;
	}
	inline T operator [] (const Vector<T> &u, const Vector<T> &v, const Vector<T> &w) {
		T r = 0;
		if (u.d != 3 || v.d != 3 || w.d != 3) return r;
		r = (*u.data)*(v.data[1] * w.data[2] - w.data[1] * v.data[2])
			- (*v.data)*(u.data[1] * w.data[2] - w.data[1] * u.data[2])
			+ (*w.data)*(u.data[1] * v.data[2] - v.data[1] * u.data[2]);
		return r;
	}*/

	friend Vector<T> operator * (matrix<T> &M, Vector<T> &v) {
		Vector<T> r(v.d);
		if (!M.square() || M.height() != v.d) return Vector();
		for (int i = 0; i < v.d; i++) {
			for (int j = 0; j < v.d; j++) {
				r.data[i] += M[i][j] * v.data[j];
			}
		}
		return r;
	}

};

#endif

class complex {
#define _i complex(0,1)
	friend inline complex timesi(const complex &a) {
		return complex(-a.ima, a.rel);
	}
	friend inline complex times_i(const complex &a) {
		return complex(a.ima, -a.rel);
	}
public:
	double rel, ima;
	complex() : rel(0), ima(0) {}
	complex(double real, double imag) : rel(real), ima(imag) {}
	complex(const complex& a) : rel(a.rel), ima(a.ima) {}
	complex(const int &a) : rel(a), ima(0) {}
	complex(const double &a) : rel(a), ima(0) {}

	friend inline complex operator + (const complex &a, const complex &b) {
		return complex(a.rel + b.rel, a.ima + b.ima);
	}
	friend inline void operator += (complex &a, const complex &b) {
		a.rel += b.rel, a.ima += b.ima;
	}
	friend inline complex operator + (const complex &a, const double &b) {
		return complex(a.rel + b, a.ima);
	}
	friend inline complex operator + (const double &a, const complex &b) {
		return complex(b.rel + a, b.ima);
	}
	friend inline void operator += (complex &a, const double &b) {
		a.rel += b;
	}

	friend inline complex operator - (const complex &a) {
		return complex(-a.rel, -a.ima);
	}
	friend inline complex operator - (const complex &a, const complex &b) {
		return complex(a.rel - b.rel, a.ima - b.ima);
	}
	friend inline complex operator - (const complex &a, const double &b) {
		return complex(a.rel - b, a.ima);
	}
	friend inline complex operator - (const double &a, const complex &b) {
		return complex(a - b.rel, -b.ima);
	}
	friend inline void operator -= (complex &a, const double &b) {
		a.rel -= b;
	}
	friend inline void operator -= (complex &a, const complex &b) {
		a.rel -= b.rel, a.ima -= b.ima;
	}

	friend inline complex operator * (const complex &a, const complex &b) {
		return complex(a.rel*b.rel - a.ima*b.ima, a.rel*b.ima + a.ima*b.rel);
	}
	friend inline void operator *= (complex &a, const complex &b) {
		a = a * b;
	}
	friend inline complex operator * (const complex &a, const double &b) {
		return complex(a.rel * b, a.ima * b);
	}
	friend inline complex operator * (const double &a, const complex &b) {
		return complex(b.rel * a, b.ima * a);
	}
	friend inline void operator *= (complex &a, const double &b) {
		a.rel *= b, a.ima *= b;
	}

	friend inline complex operator / (const complex &a, const double &b) {
		return complex(a.rel / b, a.ima / b);
	}
	friend inline complex operator / (const complex &a, const complex &b) {
		return complex(a.rel*b.rel + a.ima*b.ima, a.ima*b.rel - a.rel*b.ima) / (b.rel*b.rel + b.ima*b.ima);
	}
	friend inline void operator /= (complex &a, const complex &b) {
		a = a / b;
	}
	friend inline complex operator / (const double &a, const complex &b) {
		return complex(b.rel, -b.ima) / (a*(b.rel*b.rel + b.ima*b.ima));
	}
	friend inline void operator /= (complex &a, const double &b) {
		a.rel /= b, a.ima /= b;
	}

	friend inline double abs(complex a) {
		return sqrt(a.rel*a.rel + a.ima*a.ima)/**((a.rel < 0) ^ (a.ima < 0) ? -1 : 1)*/;
	}
	friend inline double arg(complex a) {
		return atan2(a.ima, a.rel);
		/*if (a.rel > 0) return atan(a.ima / a.rel);
		if (a.rel < 0) {
			if (a.ima >= 0) return atan(a.ima / a.rel) + _const_pi;
			return atan(a.ima / a.rel) - _const_pi;
		}
		if (a.ima > 0) return 1.57079632679489661923;
		if (a.ima < 0) return -1.57079632679489661923;
		return 0;*/
	}

	friend inline complex sqrt(const complex &a) {
		double m = abs(a);
		if (a.ima != 0) return complex(sqrt((m + a.rel) / 2), sgn(a.ima)*sqrt((m - a.rel) / 2));
		else if (a.rel >= 0) return complex(sqrt(a.rel));
		else return complex(0, sqrt(-a.rel));
	}
	friend inline complex sqrti(double a) {
		return (a >= 0) ? complex(sqrt(a)) : complex(0, sqrt(a));
	}
	friend inline complex square(const complex &a) {
		return complex(a.rel*a.rel - a.ima*a.ima, 2 * a.rel*a.ima);
	}
	friend inline complex cube(const complex &a) {
		return complex(a.rel*a.rel*a.rel - 3 * a.rel*a.ima*a.ima, 3 * a.rel*a.rel*a.ima - a.ima*a.ima*a.ima);
	}
	friend inline complex sqr4(const complex &a) {
		return complex(pow(a.rel, 4) - 6 * pow(a.rel*a.ima, 2) + pow(a.ima, 4), 4 * (pow(a.rel, 3)*a.ima - a.rel*pow(a.ima, 3)));
	}
	friend inline complex sqr5(const complex &a) {
		return complex(pow(a.rel, 5) - 10 * a.rel*a.rel*a.rel*a.ima*a.ima + 5 * a.rel*pow(a.ima, 4),
			5 * pow(a.rel, 4)*a.ima - 10 * a.rel*a.rel*a.ima*a.ima*a.ima + pow(a.ima, 5));
	}
	friend inline complex sin(const complex &a) {
		return complex(sin(a.rel)*cosh(a.ima), cos(a.rel)*sinh(a.ima));
	}
	friend inline complex cos(const complex &a) {
		return complex(cos(a.rel)*cosh(a.ima), -sin(a.rel)*sinh(a.ima));
	}
	friend inline complex tan(const complex &a) {
		return sin(a) / cos(a);
	}
	friend inline complex exp(const complex &a) {
		return exp(a.rel)*complex(cos(a.ima), sin(a.ima));
	}
	friend inline complex sinh(const complex &a) {
		return complex(sinh(a.rel)*cos(a.ima), cosh(a.rel)*sin(a.ima));
	}
	friend inline complex cosh(const complex &a) {
		return complex(cosh(a.rel)*cos(a.ima), sinh(a.rel)*sin(a.ima));
	}
	friend inline complex tanh(const complex &a) {
		return sinh(a) / cosh(a);
	}
	friend inline complex log(const complex &a) {
		return complex(log(a.rel*a.rel + a.ima*a.ima) / 2, arg(a));
	}
	friend inline complex log(const complex &a, const complex &b) {
		return log(a) / log(b);		// Could be simplified. 
	}
	friend inline complex log(const complex &a, const double &b) {
		return log(b) / log(a);
	}
	friend inline complex log(const double &a, const complex &b) {
		return log(b) / log(a);
	}
	friend inline complex log2(const complex &a) {
		return log(a) / 0.6931471805599453;
	}
	friend inline complex log10(const complex &a) {
		return log(a) / 2.302585092994046;
	}
	friend inline complex pow(const complex &a, const complex &b) {
		double sqs = a.rel*a.rel + a.ima*a.ima;
		double Arg = atan2(a.ima, a.rel);
		double inr = b.rel*Arg + b.ima*log(sqs) / 2;
		double cof = pow(sqs, b.rel / 2)*exp(-b.ima*Arg);
		return complex(cof*cos(inr), cof*sin(inr));
		//return exp(b*log(a));
	}
	friend inline complex pow(const complex &a, const double &b) {
		//return exp(b*log(a));
		double m = b * arg(a);
		return pow(a.rel*a.rel + a.ima*a.ima, b / 2) * complex(cos(m), sin(m));
	}
	friend inline complex pow(const double &a, const complex &b) {
		double inr = b.ima*log(a*a) / 2;
		double cof = pow(a*a, b.rel / 2);
		return complex(cof*cos(inr), cof*sin(inr));
		//return exp(b*log(a));
	}
	friend inline complex asin(complex a) {
		return times_i(log(complex(-a.ima, a.rel) + sqrt(1 - square(a))));
	}
	friend inline complex acos(complex a) {
		return times_i(log(a + timesi(sqrt(1 - square(a)))));
	}
	friend inline complex atan(complex a) {
		a = log(complex(-a.rel, 1 - a.ima) / complex(a.rel, a.ima + 1)) / 2;
		return complex(a.ima, -a.rel);
	}
	friend inline complex acoti(complex a) {
		a = log(complex(a.rel, 1 + a.ima) / complex(a.rel, a.ima - 1)) / 2;
		return complex(a.ima, -a.rel);
	}
	friend inline complex asinh(complex a) {
		return log(a + sqrt(square(a) + 1));
	}
	friend inline complex acosh(complex a) {
		return log(a + sqrt(square(a) - 1));
	}
	friend inline complex atanh(complex a) {
		return 0.5*log((1 + a) / (1 - a));
	}

	friend inline complex tgamma(complex z) {
		complex c = 1 / z;
		c *= pow(100, z);
		for (double n = 1; n < 100; n++) {
			//c *= pow(1 + 1 / n, z) / (1 + z / n);
			c /= (1 + z / n);
		}
		return c;
	}
	friend inline complex Euler(complex q) {
		complex c = 1;
		for (double k = 1; k < 100; k++) {
			c *= (1 - pow(q, k));
		}
		return c;
	}
	friend inline complex erf(complex a) {
		/* Taylor Series, polynomial function don't work very well for large number. */
		complex c = 0;
		for (double Pn = 1, n = 1; n < 100; n++, Pn *= n) {		// Maximum 170! 
			c -= pow(a, 2 * n) / (Pn*(2 * n + 1));
			n++, Pn *= n;
			c += pow(a, 2 * n) / (Pn*(2 * n + 1));
		}
		c += 1;
		return 2 * a / sqrt(_const_pi)*c;
	}

	friend inline complex Mandelbrot(complex c) {
		complex z = 0;
		for (int i = 0; i < 100; i++) {
			z = square(z) + c;
			//if (baddouble(z.rel) || baddouble(z.ima)) break;
		}
		return z;
	}

};
