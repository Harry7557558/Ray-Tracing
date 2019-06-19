#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
using namespace std;

typedef unsigned char byte;

#define pause system("pause");

#ifndef PI
#define PI 3.1415926535897932384626433832795029
#endif

#ifndef _INC_BitMap

#define _INC_BitMap

enum WebSafeColour {
	AliceBlue = 0xF0F8FF, AntiqueWhite = 0xFAEBD7, Aqua = 0x00FFFF, Aquamarine = 0x7FFFD4, Azure = 0xF0FFFF, Beige = 0xF5F5DC, Bisque = 0xFFE4C4, Black = 0x000000,
	BlanchedAlmond = 0xFFEBCD, Blue = 0x0000FF, BlueViolet = 0x8A2BE2, Brown = 0xA52A2A, BurlyWood = 0xDEB887, CadetBlue = 0x5F9EA0, Chartreuse = 0x7FFF00, Chocolate = 0xD2691E,
	Coral = 0xFF7F50, CornflowerBlue = 0x6495ED, Cornsilk = 0xFFF8DC, Crimson = 0xDC143C, Cyan = 0x00FFFF, DarkBlue = 0x00008B, DarkCyan = 0x008B8B, DarkGoldenRod = 0xB8860B,
	DarkGray = 0xA9A9A9, DarkGrey = 0xA9A9A9, DarkGreen = 0x006400, DarkKhaki = 0xBDB76B, DarkMagenta = 0x8B008B, DarkOliveGreen = 0x556B2F, DarkOrange = 0xFF8C00, DarkOrchid = 0x9932CC,
	DarkRed = 0x8B0000, DarkSalmon = 0xE9967A, DarkSeaGreen = 0x8FBC8F, DarkSlateBlue = 0x483D8B, DarkSlateGray = 0x2F4F4F, DarkSlateGrey = 0x2F4F4F, DarkTurquoise = 0x00CED1, DarkViolet = 0x9400D3,
	DeepPink = 0xFF1493, DeepSkyBlue = 0x00BFFF, DimGray = 0x696969, DimGrey = 0x696969, DodgerBlue = 0x1E90FF, FireBrick = 0xB22222, FloralWhite = 0xFFFAF0, ForestGreen = 0x228B22,
	Fuchsia = 0xFF00FF, Gainsboro = 0xDCDCDC, GhostWhite = 0xF8F8FF, Gold = 0xFFD700, GoldenRod = 0xDAA520, Gray = 0x808080, Grey = 0x808080, Green = 0x008000,
	GreenYellow = 0xADFF2F, HoneyDew = 0xF0FFF0, HotPink = 0xFF69B4, IndianRed = 0xCD5C5C, Indigo = 0x4B0082, Ivory = 0xFFFFF0, Khaki = 0xF0E68C, Lavender = 0xE6E6FA,
	LavenderBlush = 0xFFF0F5, LawnGreen = 0x7CFC00, LemonChiffon = 0xFFFACD, LightBlue = 0xADD8E6, LightCoral = 0xF08080, LightCyan = 0xE0FFFF, LightGoldenRodYellow = 0xFAFAD2, LightGray = 0xD3D3D3,
	LightGrey = 0xD3D3D3, LightGreen = 0x90EE90, LightPink = 0xFFB6C1, LightSalmon = 0xFFA07A, LightSeaGreen = 0x20B2AA, LightSkyBlue = 0x87CEFA, LightSlateGray = 0x778899, LightSlateGrey = 0x778899,
	LightSteelBlue = 0xB0C4DE, LightYellow = 0xFFFFE0, Lime = 0x00FF00, LimeGreen = 0x32CD32, Linen = 0xFAF0E6, Magenta = 0xFF00FF, Maroon = 0x800000, MediumAquaMarine = 0x66CDAA,
	MediumBlue = 0x0000CD, MediumOrchid = 0xBA55D3, MediumPurple = 0x9370DB, MediumSeaGreen = 0x3CB371, MediumSlateBlue = 0x7B68EE, MediumSpringGreen = 0x00FA9A, MediumTurquoise = 0x48D1CC, MediumVioletRed = 0xC71585,
	MidnightBlue = 0x191970, MintCream = 0xF5FFFA, MistyRose = 0xFFE4E1, Moccasin = 0xFFE4B5, NavajoWhite = 0xFFDEAD, Navy = 0x000080, OldLace = 0xFDF5E6, Olive = 0x808000,
	OliveDrab = 0x6B8E23, Orange = 0xFFA500, OrangeRed = 0xFF4500, Orchid = 0xDA70D6, PaleGoldenRod = 0xEEE8AA, PaleGreen = 0x98FB98, PaleTurquoise = 0xAFEEEE, PaleVioletRed = 0xDB7093,
	PapayaWhip = 0xFFEFD5, PeachPuff = 0xFFDAB9, Peru = 0xCD853F, Pink = 0xFFC0CB, Plum = 0xDDA0DD, PowderBlue = 0xB0E0E6, Purple = 0x800080, RebeccaPurple = 0x663399,
	Red = 0xFF0000, RosyBrown = 0xBC8F8F, RoyalBlue = 0x4169E1, SaddleBrown = 0x8B4513, Salmon = 0xFA8072, SandyBrown = 0xF4A460, SeaGreen = 0x2E8B57, SeaShell = 0xFFF5EE,
	Sienna = 0xA0522D, Silver = 0xC0C0C0, SkyBlue = 0x87CEEB, SlateBlue = 0x6A5ACD, SlateGray = 0x708090, SlateGrey = 0x708090, Snow = 0xFFFAFA, SpringGreen = 0x00FF7F,
	SteelBlue = 0x4682B4, Tan = 0xD2B48C, Teal = 0x008080, Thistle = 0xD8BFD8, Tomato = 0xFF6347, Turquoise = 0x40E0D0, Violet = 0xEE82EE, Wheat = 0xF5DEB3,
	White = 0xFFFFFF, WhiteSmoke = 0xF5F5F5, Yellow = 0xFFFF00, YellowGreen = 0x9ACD32
};

class pixel {
public:
	byte b, g, r;
	friend class bitmap;
	pixel()/* : r(0), g(0), b(0) */ {}
	pixel(const byte &R, const byte &G, const byte &B) :r(R), g(G), b(B) {}
	pixel(const pixel& a) :r(a.r), g(a.g), b(a.b) {}
	pixel(unsigned hex) {
		b = hex, hex >>= 8; g = hex, hex >>= 8; r = hex;
	}
	pixel(const WebSafeColour &color) {
		b = unsigned(color), g = unsigned(color) >> 8, r = unsigned(color) >> 16;
	}
	void inverse() {
		b = ~b, g = ~g, r = ~r;
	}
	friend ostream& operator << (ostream& os, pixel n) {
		os << "#" << hex << uppercase << setw(2) << setfill('0') << (unsigned)n.r << setw(2) << setfill('0') << (unsigned)n.g << setw(2) << setfill('0') << (unsigned)n.b;
		return os;
	}
};
inline pixel rgb(const byte &r, const byte &g, const byte &b) {
	pixel k(r, g, b); return k;
}
inline pixel rgb(const byte &l) {
	pixel k(l, l, l); return k;
}
inline pixel drgb(const double &r, const double &g, const double &b) {
	pixel k(r * 255, g * 255, b * 255); return k;
}
inline pixel color(const WebSafeColour &col) {
	return pixel(byte(unsigned(col) >> 16), byte(unsigned(col) >> 8), byte(unsigned(col)));
}

inline double Myt_fromHSL_HueToRGB(double v1, double v2, double vH) {
	vH -= floor(vH);
	if ((6 * vH) < 1) return (v1 + (v2 - v1) * 6 * vH);
	if ((2 * vH) < 1) return v2;
	if ((3 * vH) < 2) return (v1 + (v2 - v1) * ((2.0 / 3) - vH) * 6);
	return v1;
}
inline pixel fromHSL(double h, double s, double l) {
	// 0 ≤ h < 1, 0 ≤ s,l ≤ 1
	pixel C;
	if (s == 0)
	{
		C.r = C.g = C.b = l * 255;
		return C;
	}
	double v1, v2;
	v2 = (l < 0.5) ? (l * (1 + s)) : ((l + s) - (l * s));
	v1 = 2 * l - v2;
	C.r = (unsigned char)(255 * Myt_fromHSL_HueToRGB(v1, v2, h + (1.0f / 3)));
	C.g = (unsigned char)(255 * Myt_fromHSL_HueToRGB(v1, v2, h));
	C.b = (unsigned char)(255 * Myt_fromHSL_HueToRGB(v1, v2, h - (1.0f / 3)));
	return C;
}


class bitmap {

protected:
	unsigned w, h;	// width and height
	unsigned bpp;		// bits per pixel
	unsigned xpm, ypm;	// X and Y pixels per meter
	pixel *data; 	// data

	friend class GIF;

public:
	bitmap() : w(0), h(0), bpp(0), xpm(0), ypm(0), data(0) {}
	bitmap(unsigned width, unsigned height) : w(width), h(height) {
		bpp = 24;
		xpm = ypm = 0;
		data = new pixel[w*h];
		for (int i = w * h; i > 0; i--) data->r = data->g = data->b = 255, data++;
		data -= w * h;
	}
	bitmap(unsigned width, unsigned height, pixel color) : w(width), h(height) {
		bpp = 24;
		xpm = ypm = 0;
		int n = w * h;
		data = new pixel[w*h];
		for (int i = 0; i < n; i++) *data = color, data++;
		data -= n;
	}
	bitmap(string FilePath) {
		w = h = 0, data = 0;
		ifstream _Image(FilePath, ios_base::in | ios_base::binary);
		if (_Image.fail()) return;
		char f[4];
		_Image.read(f, 2);	// Signature, BM
		if (f[0] != 'B' || f[1] != 'M') return;
		unsigned Fsize;
		_Image.read(f, 4);	// File Size
		Fsize = *((unsigned*)f);
		_Image.read(f, 4);	// reserved
		_Image.read(f, 4);	// Full header size
		if (*((unsigned*)f) != 0x36) return;	// ⚠
		_Image.read(f, 4);	// DIB header size
		if (*((unsigned*)f) != 0x28) return;	// ⚠
		_Image.read(f, 4);	// width
		w = *((unsigned*)f);
		_Image.read(f, 4);	// height
		h = *((unsigned*)f);
		_Image.read(f, 2);	// planes, always 1
		if (*((unsigned short*)f) != 1) return;
		_Image.read(f, 2);	// bits per pixel, 8 bit RGB is 0x18
		bpp = *((unsigned short*)f);
		if (bpp != 24 && bpp != 32) return;	// ⚠
		_Image.read(f, 4);	// compression
		if (*((unsigned*)f) != 0) return;	// ⚠
		_Image.read(f, 4);	// Image size ⚠
		_Image.read(f, 4);	// X pixels per meter
		xpm = *((unsigned*)f);
		_Image.read(f, 4);	// Y pixels per meter
		ypm = *((unsigned*)f);
		_Image.read(f, 4);
		if (*((unsigned*)f) != 0) return;	// ⚠
		_Image.read(f, 4);
		if (*((unsigned*)f) != 0) return;	// ⚠

		_Image.seekg(0, ios::end);
		Fsize -= 54;	// Given image size
		unsigned Rsize = _Image.tellg(); Rsize -= 54;	// Measured image size
		int padding = (4 - ((w * bpp / 8) % 4)) % 4;
		unsigned Csize = w * h*bpp / 8 + h * padding;	// Theoretical image size, follow the standard
		bool pad = 1;
		if (Fsize != Rsize) return;
		if (Fsize == Csize);
		else if (Fsize - Csize == 2) Fsize -= 2;	// GAP, optional
		else if (Fsize - Csize == h * padding) pad = 0;	// Lose padding
		else if (Fsize - Csize - 2 == h * padding) pad = 0;	// GAP + No Padding
		else return;

		_Image.seekg(54);
		data = new pixel[w*h];
		switch (bpp) {
		case 24: {
			for (unsigned i = 0; i < h; i++) {
				for (unsigned j = 0; j < w; j++) {
					_Image.read(f, 3);
					data->b = f[0], data->g = f[1], data->r = f[2];
					data++;
				}
				if (pad) _Image.read(f, padding);
			}
			break;
		}
		case 32: {
			for (unsigned i = 0; i < h; i++) {
				for (unsigned j = 0; j < w; j++) {
					_Image.read(f, 4);
					data->b = f[0], data->g = f[1], data->r = f[2];
					data++;
				}
			}
			bpp = 24;
			break;
		}
		default: {
			delete data;
			return;
		}
		}
		data -= w * h;
		_Image.close();
	}
	~bitmap() {
		if (data != 0) delete data;
		w = h = 0;
		data = 0;
	}
	bitmap(const bitmap &a) {
		w = a.w, h = a.h;
		bpp = a.bpp, xpm = a.xpm, ypm = a.ypm;
		data = new pixel[w*h];
		const pixel* p = a.data;
		for (int i = w * h; i > 0; i--) {
			*data = *p; data++, p++;
		}
		data -= w * h;
	}
	bitmap& operator = (const bitmap &a) {
		if (this->data != 0) delete this->data;
		w = a.w, h = a.h;
		bpp = a.bpp, xpm = a.xpm, ypm = a.ypm;
		data = new pixel[w*h];
		const pixel* p = a.data;
		for (int i = w * h; i > 0; i--) {
			*data = *p; data++, p++;
		}
		data -= w * h;
		return *this;
	}
	inline bool fail() {
		return (data == 0 || w == 0 || h == 0 || bpp == 0);
	}
	bool out(string FilePath) {
		if (data == 0) return 0;
		if (w == 0 || h == 0) return 0;
		//bpp = 0x18;
		int padding = 4 - (w * bpp / 8) % 4;
		if (w*bpp % 8 != 0) padding--;
		padding %= 4;
		ofstream _Image(FilePath, ios_base::out | ios_base::binary);
		if (_Image.fail()) {
			cout << "\aCreate file \"" << FilePath << "\" failed. \n";
			cout << "Press any key to retry. \n";
			pause;
			_Image.open(FilePath, ios_base::out | ios_base::binary);
			if (_Image.fail()) {
				cout << "\aCreate file \"" << FilePath << "\" failed. \n";
				cout << "Press any key to retry. \n";
				pause;
				_Image.open(FilePath, ios_base::out | ios_base::binary);
				if (_Image.fail()) {
					cout << "\aCreate file \"" << FilePath << "\" failed. \n";
					cout << "Press any key to retry. \n";
					pause;
					_Image.open(FilePath, ios_base::out | ios_base::binary);
					if (_Image.fail()) {
						cout << "\aCreate file \"" << FilePath << "\" failed. \n";
						cout << "Press Enter to exit . . . ";
						getchar();
						return false;
					}
				}
			}
		}
		_Image << "BM";	// Signature, BM
		unsigned size = w * h * bpp / 8 + h * padding + 54;
		unsigned os;
		_Image.write((char*)&size, 4);	// File Size
		os = 0;
		_Image.write((char*)&os, 4);	// reserved
		os = bpp == 1 ? 0x3E : 0x36;
		_Image.write((char*)&os, 4);	// Full header size, usually 54 byte
		os = 0x28;
		_Image.write((char*)&os, 4);	// DIB header size, 40 byte
		_Image.write((char*)&w, 4);	// width
		_Image.write((char*)&h, 4);	// height
		os = 1;
		_Image.write((char*)&os, 2);	// planes, always 1
		_Image.write((char*)&bpp, 2);	// bits per pixel ⚠
		os = 0;
		_Image.write((char*)&os, 4);	// compression
		size -= 54;
		_Image.write((char*)&size, 4);	// Image Size
		_Image.write((char*)&xpm, 4);	// X pixels per meter
		_Image.write((char*)&ypm, 4);	// Y pixels per meter
		_Image.write((char*)&os, 4);
		_Image.write((char*)&os, 4);
		size = 0;
		pixel *b = data;
		if (bpp == 0x18) {
			for (unsigned i = 0; i < h; i++) {
				for (unsigned j = 0; j < w; j++) {
					_Image.write((char*)b, 3);
					//if (b->r != 0) cout << "[" << i << "][" << j << "]" << int(b->r) << " ";
					b++;
				}
				_Image.write((char*)&size, padding);	// Pad row size to a multiple of 4
			}
		}
		else if (bpp == 1) {
			_Image.write((char*)&os, 4);
			os = 0xFFFFFF;
			_Image.write((char*)&os, 4);
			byte s;
			for (int i = 0; i < h; i++) {
				os = 0; s = 0;
				for (int j = 0; j < w; j++) {
					s <<= 1;
					if (b->r != 0) s |= 1;
					if (++os == 8) {
						_Image.write((char*)&s, 1);
						os = 0;
						s = 0;
					}
					b++;
				}
				if (os != 0) _Image.write((char*)&s, 1);
				os = 0;
				s = 0;
				_Image.write((char*)&os, padding);	// Pad row size to a multiple of 4
			}
		}
		b = 0;
		_Image.close();
		return 1;
	}

	inline constexpr unsigned height() const { return h; }
	inline constexpr unsigned width() const { return w; }
	inline constexpr unsigned size() const { return h * w; }
	inline pixel* operator [] (int k) {
		return &data[k*w];
	}
	void clear() {
		pixel *b = data;
		for (int i = w * h; i > 0; i--) {
			b->b = b->g = b->r = 0;
			b++;
		}
	}
	void clear(unsigned width, unsigned height) {
		delete data;
		data = new pixel[width*height];
		w = width, h = height;
	}
	void clear(pixel a) {
		pixel *b = data;
		for (int i = w * h; i > 0; i--) {
			*b = a;
			b++;
		}
	}
	void clear(unsigned width, unsigned height, pixel a) {
		delete[] data;
		data = new pixel[width*height];
		w = width, h = height;
		pixel *b = data;
		for (int i = w * h; i > 0; i--) {
			*b = a; b++;
		}
	}
	bitmap subimg(int x0, int y0, unsigned x, unsigned y) const {
		x += x0, y += y0;
		if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
		if (x > w) x = w; if (y > h) y = h;
		bitmap r;
		r.w = x - x0, r.h = y - y0;
		r.bpp = this->bpp, r.xpm = this->xpm, r.ypm = this->ypm;
		r.data = new pixel[r.w*r.h];
		const pixel* p;
		for (int i = y0; i < y; i++) {
			p = this->data + i * w + x0;
			for (int j = x0; j < x; j++) {
				*r.data = *p; r.data++, p++;
			}
		}
		r.data -= r.w*r.h;
		return r;
	}
	void setbpp(unsigned n) {
		if (n == 24) bpp = 24;
		else if (n == 32) bpp = 32;
		else if (n == 1) this->binarize();
		else;
	}

	// Drawing Patterns
	inline void dot(int x, int y, pixel fill) {
		if (x >= 0 && x < w && y >= 0 && y < h) data[y * w + x] = fill;
	}
	inline void dot(int x, int y, WebSafeColour fill) {
		if (x >= 0 && x < w && y >= 0 && y < h) data[y * w + x] = fill;
	}
	inline void dot(int x, int y, pixel fill, double alpha) {
		if (x >= 0 && x < w && y >= 0 && y < w) {
			pixel *p = data + y * w + x;
			*p = rgb(fill.r + (p->r - double(fill.r))*alpha, fill.g + (p->g - double(fill.g))*alpha, fill.b + (p->b*double(fill.b))*alpha);
		}
	}
	void rect(int x, int y, int width, int height, pixel fill) {
		if (width == 0 || height == 0) return;
		if (width < 0) x += width, width = -width;
		if (height < 0) y += height, height = -height;
		if (x < 0) width += x, x = 0;
		if (y < 0) height += y, y = 0;
		if (width < 0 || height < 0) return;
		if (x > w || y > h) return;
		if (x + width > w) width -= x + width - w;
		if (y + height > h) height -= y + height - h;
		pixel *p = data + y * w + x;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				*p = fill; p++;
			}
			p += w - width;
		}
	}
	void rect(int x, int y, int width, int height, pixel fill, double alpha) {
		if (width == 0 || height == 0) return;
		if (width < 0) x += width, width = -width;
		if (height < 0) y += height, height = -height;
		if (x < 0) width += x, x = 0;
		if (y < 0) height += y, y = 0;
		if (width < 0 || height < 0) return;
		if (x > w || y > h) return;
		if (x + width > w) width -= x + width - w;
		if (y + height > h) height -= y + height - h;
		pixel *p = data + y * w + x;
		alpha = 1 - alpha;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				*p = rgb(fill.r + (p->r - fill.r)*alpha, fill.g + (p->g - fill.g)*alpha, fill.b + (p->b*fill.b)*alpha);
				p++;
			}
			p += w - width;
		}
	}
	void circle(int cx, int cy, int r, pixel fill) {
		if (r < 0) r = -r;
		if (r == 0) return;
		cx -= r, cy -= r;
		int d = 2 * r;
		pixel *p;
		double rd;
		for (int i = 0; i < d; i++) {
			rd = sqrt(r*r - (r - i)*(r - i)) * 2;
			p = data + (cy + i) * w + int(r - rd / 2) + cx;
			for (int j = 0; j < rd; j++) {
				*p = fill, p++;
			}
		}
	}
	void line(int x1, int y1, int x2, int y2, double width, pixel fill) {
		// https://zhuanlan.zhihu.com/p/30553006
		/*if (x2 == x1) {
			rect(x1, y1 - width / 2, width, y2 - y1, fill);
			return;
		}*/
		double slope, intercept;
		pixel *p; int d;
		if (abs(y2 - y1) < abs(x2 - x1)) {
			if (x1 > x2) swap(x1, x2), swap(y1, y2);
			if (x2 < 0 || (y1 < 0 && y2 < 0)) return;
			if (x1 > int(w) || (y1 > int(h) && y2 > int(h))) return;
			slope = double(y2 - y1) / (x2 - x1);
			intercept = y1 - x1 * slope;
			//int swidth = double(width / slope) * sqrt(slope*slope + 1);
			//int hwidth = width * sqrt(slope*slope + 1);
			if (x1 < 0) x1 = 0, y1 = intercept;
			if (x2 > w) x2 = w, y2 = w * slope + intercept;
			if (x1*slope + intercept < 0) {
				x1 = -intercept / slope;
				y1 = 0;
			}
			else if (x1*slope + intercept > h) {
				x1 = (h - intercept) / slope;
				y1 = h;
			}
			if (x2*slope + intercept < 0) {
				x2 = -intercept / slope;
				y2 = 0;
			}
			else if (x2*slope + intercept > h) {
				x2 = (h - intercept) / slope;
				y1 = h;
			}
			//y1 += hwidth / 2, y2 += hwidth / 2, intercept += hwidth / 2;
			for (int i = x1; i < x2; i++) {
				d = i * slope + intercept;
				if (d >= 0 && d < h) {
					p = data + d * w + i;
					//*p = rgb(fill.r + (p->r - fill.r) / 2, fill.g + (p->g - fill.g) / 2, fill.b + (p->b*fill.b) / 2), p++;
					//for (int j = 2; j < swidth; j++) {
					*p = fill, p++;
					//}
					//*p = rgb(fill.r + (p->r - fill.r) / 2, fill.g + (p->g - fill.g) / 2, fill.b + (p->b*fill.b) / 2), p++;
				}
			}
		}
		else {
			if (y1 > y2) swap(x1, x2), swap(y1, y2);
			if (y2 < 0 || (x1 < 0 && x2 < 0)) return;
			if (y1 > int(h) || (x1 > int(w) && x2 > int(w))) return;
			slope = double(x2 - x1) / (y2 - y1);
			intercept = x1 - y1 * slope;
			if (y1 < 0) y1 = 0, x1 = intercept;
			if (y2 > h) y2 = h, x2 = h * slope + intercept;
			if (y1*slope + intercept < 0) {
				y1 = -intercept / slope;
				x1 = 0;
			}
			else if (y1*slope + intercept > w) {
				y1 = (w - intercept) / slope;
				x1 = w;
			}
			if (y2*slope + intercept < 0) {
				y2 = -intercept / slope;
				x2 = 0;
			}
			else if (y2*slope + intercept > w) {
				y2 = (w - intercept) / slope;
				x1 = w;
			}
			for (int i = y1; i < y2; i++) {
				d = i * slope + intercept;
				if (d >= 0 && d < w) {
					p = data + i * w + d;
					*p = fill, p++;
				}
			}
		}
	}

	// Filters
	void binarize() {
		const unsigned size = w * h;
		for (int i = 0; i < size; i++) {
			//if ((unsigned(data->r) * 19595 + unsigned(data->g) * 38469 + unsigned(data->b) * 7472) > 0x1800000) data->r = data->g = data->b = 1;
			if (0.299*data->r + 0.587*data->g + 0.114*data->b > 128) data->r = data->g = data->b = 255;
			else data->r = data->g = data->b = 0;
			data++;
		}
		data -= size;
		bpp = 1;
	}
	// Gaussian Blur: May be low efficiency
	void GaussianBlur(int r) {
		if (r < 1) return;
		float rs, gs, bs;
		pixel* p = this->data, *q;
		bitmap M(*this);
		int x0, y0, x1, y1;
		unsigned rt = r; r = ceil(r*2.57);
		float *ts = new float[4 * r*r];
		for (int i = -r; i < r; i++) {
			for (int j = -r; j < r; j++) {
				*ts = exp((i*i + j * j) / (-2.0 * rt*rt)) / (2 * PI*rt*rt);
				ts++;
			}
		}
		ts -= 4 * r * r;

		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				rs = gs = bs = 0;
				x0 = j - r, y0 = i - r, x1 = j + r, y1 = i + r;
				for (int m = y0; m < y1; m++) {
					q = data + m * w + x0;
					for (int n = x0; n < x1; n++) {
						if (m > 0 && n > 0 && m < h && n < w) {
							rs += (*ts)*q->r, gs += (*ts)*q->g, bs += (*ts)*q->b;
						}
						q++, ts++;
					}
				}
				ts -= 4 * r * r;

				p->r = rs, p->g = gs, p->b = bs;
				p++;
			}
		}
		p = 0;
	}
};

#endif

