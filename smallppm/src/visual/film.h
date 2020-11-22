#pragma once

#include "math/linagl.h"
#include "math/filter.h"
#include "math/math_utils.h"
#include "system/threading.h"
#include <mutex>
#include <vector>
#include <algorithm>
#include <fstream>

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
	Film(int w, int h, Filter* pFilter = nullptr) : resX(w), resY(h) {
		aspect = (real)(resX) / (real)(resY);
		//filter = std::unique_ptr<BoxFilter>(new BoxFilter());
		pixelBuffer.resize(resX * resY);
		imageBuffer.resize(resX * resY);
		bufferLocks.resize(resX * resY);
		if (pFilter == nullptr) {
			filter.reset(new BoxFilter());
		}
		else {
			filter.reset(pFilter);
		}

	}
	real Area() const {
		//return area;
		return width * height;
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
				std::lock_guard<Spinlock> lock(bufferLocks[pixelIndex]);
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
		//imageBuffer = image;
		for (int i = 0; i < resY; ++i) {
			for (int j = 0; j < resX; ++j) {
				int index = i * resX + j;
				Pixel& pixel = pixelBuffer[index];
				pixel.color = image[index];
				pixel.weight = 1;
			}
		}
	}

	void SetFileName(const std::string &pFileName) {
		filename = pFileName;
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
				Pixel &pixel = pixelBuffer[index];
				imageBuffer[index] = pixel.color / pixel.weight + pixel.splat;
			}
		}
	}


private:
	std::string filename;
	std::vector<Pixel> pixelBuffer;
	std::vector<Vec3> imageBuffer;
	std::unique_ptr<Filter> filter;
	std::vector<Spinlock> bufferLocks;
};

NAMESPACE_END