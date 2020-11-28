#pragma once

#include "math/Linagl.h"
#include "math/Filter.h"
#include "math/MathUtils.h"
#include "system/Threading.h"
#include <mutex>
#include <vector>
#include <algorithm>
#include <fstream>

GY_NAMESPACE_BEGIN

struct Pixel
{
	Vec3 color;
	Vec3 splat;
	real weight;
};

class Film 
{
private:
	std::string mFilename;
	std::vector<Pixel> mPixelBuffer;
	std::vector<Vec3> mImageBuffer;
	std::unique_ptr<Filter> mpFilter;
	std::vector<Spinlock> mBufferLocks;
public:
	Film(int w, int h, Filter* pFilter = nullptr) : resX(w), resY(h) {
		aspect = (real)(resX) / (real)(resY);
		//filter = std::unique_ptr<BoxFilter>(new BoxFilter());
		mPixelBuffer.resize((int64)resX * resY);
		mImageBuffer.resize((int64)resX * resY);
		mBufferLocks.resize((int64)resX * resY);
		if (pFilter == nullptr) {
			mpFilter.reset(new BoxFilter());
		}
		else {
			mpFilter.reset(pFilter);
		}

	}
	real Area() const {
		//return area;
		return width * height;
	}
	void AddSample(real x, real y, const Vec3 &sample) {
		x -= 0.5;
		y -= 0.5;

		int minX = (int)(std::ceil(x - mpFilter->GetRadius()));
		int maxX = (int)(std::floor(x + mpFilter->GetRadius()));
		int minY = (int)(std::ceil(y - mpFilter->GetRadius()));
		int maxY = (int)(std::floor(y + mpFilter->GetRadius()));
		minX = std::max(0, minX);
		maxX = std::min(maxX, resX - 1);
		minY = std::max(0, minY);
		maxY = std::min(maxY, resY - 1);

		for (int i = minY; i <= maxY; i++) {
			for (int j = minX; j <= maxX; j++) {
				//int rowAdd = resY - 1 - i;
				//int colAdd = j;
				int pixelIndex = i * resX + j;
				std::lock_guard<Spinlock> lock(mBufferLocks[pixelIndex]);
				Pixel& pixel = mPixelBuffer[pixelIndex];
				real weight = mpFilter->Evaluate(j - x, i - y);
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
		Pixel& pixel = mPixelBuffer[pixelIndex];
		pixel.splat = pixel.splat + sample;
	}

	void SetImage(const std::vector<Vec3> &image) {
		//imageBuffer = image;
		for (int i = 0; i < resY; ++i) {
			for (int j = 0; j < resX; ++j) {
				int index = i * resX + j;
				Pixel& pixel = mPixelBuffer[index];
				pixel.color = image[index];
				pixel.weight = 1;
			}
		}
	}

	void SetFileName(const std::string &pFileName) {
		mFilename = pFileName;
	}

	void SaveImage();



public:
	int resX, resY;
	real width, height;
	real aspect;
	real area;
	Vec3 LL, LU, RL, RU;
private:

	void WriteToPixelBuffer() {
		for (int i = 0; i < resY; ++i) {
			for (int j = 0; j < resX; ++j) {
				int index = i * resX + j;
				Pixel &pixel = mPixelBuffer[index];
				mImageBuffer[index] = pixel.color / pixel.weight + pixel.splat;
			}
		}
	}

};

GY_NAMESPACE_END