#pragma once

#include "linagl.h"
#include "filter.h"
#include <algorithm>
#include <fstream>
#include "scalar.h"
#include "svpng.inc"
#include <vector>
#include "threading.h"
#include "image_io.h"
NAMESPACE_BEGIN

class Film {
protected:
	struct Pixel
	{
		Vec3 color;
		Vec3 splat;
		real weight;
	};
public:
	Film(int w, int h) : resX(w), resY(h) {
		aspect = (real)(resX) / (real)(resY);
		filter = std::unique_ptr<BoxFilter>(new BoxFilter());
	}
	real Area() const {
		return area;
	}
	void AddSample(real x, real y, const Vec3 &sample) {
		x -= 0.5;
		y -= 0.5;

		int minX = (int)(std::ceil(x - filter->GetRadius()));
		int maxX = (int)(std::floor(x + filter->GetRadius()));
		int minY = (int)(std::ceil(y - filter->GetRadius()));
		int maxY = (int)(std::floor(y + filter->GetRadius()));
		minX = std::max(0, minX);
		maxX = std::min(maxX, resX - 1);
		minY = std::max(0, minY);
		maxY = std::min(maxY, resY - 1);

		for (int i = minY; i <= maxY; i++) {
			for (int j = minX; j <= maxX; j++) {
				//int rowAdd = resY - 1 - i;
				//int colAdd = j;
				int pixelIndex = i * resX + j;
				Pixel& pixel = pixelBuffer[pixelIndex];

				real weight = filter->Evaluate(j - x, i - y);
				pixel.weight += weight;
				pixel.color = pixel.color + sample * weight;
			}
		}
	}

	void AddSplat(real x, real y, const Vec3& sample)
	{
		int X = (int)(std::floor(x));
		int Y = (int)(std::floor(y));
		X = Clamp(X, 0, resX - 1);
		Y = Clamp(Y, 0, resY - 1);

		//int rowAdd = resY - 1 - Y;
		//int colAdd = X;
		int pixelIndex = X * resX + Y;
		Pixel& pixel = pixelBuffer[pixelIndex];
		pixel.splat = pixel.splat + sample;
	}

	void SetImage(const std::vector<Vec3> &image) {
		imageBuffer = image;
	}

	void SetFileName(const std::string &pFileName) {
		filename = pFileName;
	}

	void SaveImage() {
		std::string suffix = filename.substr(filename.size() - 4);
		if (suffix == ".png") {
			ImageIO::WritePngFile(filename, imageBuffer, resX, resY);
		}
		else if (suffix == ".bmp") {
			ImageIO::WriteBmpFile(filename, imageBuffer, resX, resY);
		}
		else {

		}
	}

	/*void ConvertBmpToPng(std::string path) {

		typedef unsigned char byte;
		std::ifstream bmp(path.c_str(), std::ios::binary);
		BmpHeader header;
		char a[2];
		bmp.read(&a[0], 2);
		std::cout << a[0] << a[1] << std::endl;
		bmp.read((char*)&header, sizeof(BmpHeader));
		std::cout << header.mHeight << std::endl;
		byte *pngImage = new byte[resX * resY * 3];
		std::string output = "output.png";
		FILE* f = fopen(output.c_str(), "wb");
		for (int y = 0; y < resY; y++)
		{
			for (int x = 0; x < resX; x++)
			{

				byte bgrB[3];
				bmp.read((char*)&bgrB[0], sizeof(bgrB));
				int index = 3 * x + 3 * (resY - y - 1) * resX;
				pngImage[index] = bgrB[2];
				pngImage[index + 1] = bgrB[1];
				pngImage[index + 2] = bgrB[0];

			}
		}
		svpng(f, resX, resY, pngImage, 0);
		delete[] pngImage;
		fclose(f);
	}*/

public:
	int resX, resY;
	real width, heigh;
	real aspect;
	real area;
	Vec3 LL, LU, RL, RU;
private:
	//struct BmpHeader
	//{
	//	uint32   mFileSize;        // Size of file in bytes
	//	uint32   mReserved01;      // 2x 2 reserved bytes
	//	uint32   mDataOffset;      // Offset in bytes where data can be found (54)

	//	uint32    mHeaderSize;      // 40B
	//	uint32    mWidth;           // Width in pixels
	//	uint32    mHeight;          // Height in pixels

	//	int16  mColorPlates;     // Must be 1
	//	int16  mBitsPerPixel;    // We use 24bpp
	//	uint32   mCompression;     // We use BI_RGB ~ 0, uncompressed
	//	uint32   mImageSize;       // mWidth x mHeight x 3B
	//	uint32   mHorizRes;        // Pixels per meter (75dpi ~ 2953ppm)
	//	uint32   mVertRes;         // Pixels per meter (75dpi ~ 2953ppm)
	//	uint32   mPaletteColors;   // Not using palette - 0
	//	uint32   mImportantColors; // 0 - all are important
	//};

	//void WritePngFile() {
	//	FILE* f = fopen(filename.c_str(), "wb");
	//	typedef unsigned char byte;
	//	byte *pngImage = new byte[resX * resY * 3];
	//	for (int j = 0; j < resY; ++j) {
	//		for (int i = 0; i < resX; ++i) {
	//			int index = 3 * j * resX + 3 * i;
	//			int imageBufferIndex = j * resX + i;
	//			pngImage[index] = (byte)(toInt(imageBuffer[imageBufferIndex].x));
	//			pngImage[index + 1] = (byte)(toInt(imageBuffer[imageBufferIndex].y));
	//			pngImage[index + 2] = (byte)(toInt(imageBuffer[imageBufferIndex].z));
	//		}
	//	}
	//	svpng(f, resX, resY, pngImage, 0);
	//	delete[] pngImage;
	//	fclose(f);
	//}

	// tone mapping and gamma correction
	//int toInt(real x) {
	//	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
	//}
	//int toInt(real x) {
	//	return int(pow(Clamp(x, 0.0, 1.0), 1 / 2.2) * 255 + .5);
	//}

	void WriteToPixelBuffer() {
		for (int i = 0; i < resX; ++i) {
			for (int j = 0; j < resY; ++j) {
				int index = i * resX + j;
				Pixel &pixel = pixelBuffer[index];
				imageBuffer[index] = pixel.color / pixel.weight + pixel.splat;
			}
		}
	}

	//void WriteBmpFile()
	//{
	//	std::ofstream bmp(filename.c_str(), std::ios::binary);
	//	BmpHeader header;
	//	bmp.write("BM", 2);
	//	header.mFileSize = (uint32)(sizeof(BmpHeader) + 2) + resX * resY * 3;
	//	header.mReserved01 = 0;
	//	header.mDataOffset = (uint32)(sizeof(BmpHeader) + 2);
	//	header.mHeaderSize = 40;
	//	header.mWidth = resX;
	//	header.mHeight = resY;
	//	header.mColorPlates = 1;
	//	header.mBitsPerPixel = 24;
	//	header.mCompression = 0;
	//	header.mImageSize = resX * resY * 3;
	//	header.mHorizRes = 2953;
	//	header.mVertRes = 2953;
	//	header.mPaletteColors = 0;
	//	header.mImportantColors = 0;
	//	bmp.write((char*)&header, sizeof(header));

	//	for (int y = 0; y < resY; y++)
	//	{
	//		for (int x = 0; x < resX; x++)
	//		{
	//			const Vec &rgbF = imageBuffer[x + (resY - y - 1) * resX];
	//			typedef unsigned char byte;
	//			byte bgrB[3];
	//			bgrB[0] = byte(toInt(rgbF.z));
	//			bgrB[1] = byte(toInt(rgbF.y));
	//			bgrB[2] = byte(toInt(rgbF.x));

	//			bmp.write((char*)&bgrB, sizeof(bgrB));
	//		}
	//	}
	//}



private:
	std::string filename;
	std::vector<Pixel> pixelBuffer;
	std::vector<Vec3> imageBuffer;
	std::unique_ptr<Filter> filter;
	std::vector<Spinlock> bufferLocks;

};

NAMESPACE_END