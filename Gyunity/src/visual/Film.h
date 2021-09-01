#pragma once

#include "math/Linagl.h"
#include "math/Filter.h"
#include "math/MathUtils.h"
#include "system/Threading.h"
#include <mutex>
#include <vector>
#include <algorithm>
#include <fstream>

GYT_NAMESPACE_BEGIN

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
	int mResX, mResY;
	real mWidth, mHeight;
	real mAspectRatio;
	real mArea;
	Vec3 mLL, mLU, mRL, mRU;

public:
	Film(int w, int h, Filter* pFilter = nullptr) 
		: mResX(w), mResY(h) 
	{
		mAspectRatio = (real)(mResX) / (real)(mResY);
		//filter = std::unique_ptr<BoxFilter>(new BoxFilter());
		mPixelBuffer.resize((int64)mResX * mResY);
		mImageBuffer.resize((int64)mResX * mResY);
		mBufferLocks.resize((int64)mResX * mResY);
		if (pFilter == nullptr) {
			mpFilter.reset(new BoxFilter());
		}
		else {
			mpFilter.reset(pFilter);
		}

	}

	real Area() const 
	{
		//return area;
		return mWidth * mHeight;
	}

	void AddSample(real x, real y, const Vec3 &sample) 
	{
		x -= 0.5;
		y -= 0.5;

		int minX = (int)(std::ceil(x - mpFilter->GetRadius()));
		int maxX = (int)(std::floor(x + mpFilter->GetRadius()));
		int minY = (int)(std::ceil(y - mpFilter->GetRadius()));
		int maxY = (int)(std::floor(y + mpFilter->GetRadius()));
		minX = std::max(0, minX);
		maxX = std::min(maxX, mResX - 1);
		minY = std::max(0, minY);
		maxY = std::min(maxY, mResY - 1);

		for (int i = minY; i <= maxY; i++) {
			for (int j = minX; j <= maxX; j++) {
				//int rowAdd = resY - 1 - i;
				//int colAdd = j;
				int pixelIndex = i * mResX + j;
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
		X = Clamp(X, 0, mResX - 1);
		Y = Clamp(Y, 0, mResY - 1);

		{
			int pixelIndex = Y * mResX + X;
			std::lock_guard<Spinlock> lock(mBufferLocks[pixelIndex]);
			Pixel& pixel = mPixelBuffer[pixelIndex];
			pixel.splat = pixel.splat + sample;
		}

	}

	void SetImage(const std::vector<Vec3> &image) 
	{
		//imageBuffer = image;
		for (int i = 0; i < mResY; ++i) {
			for (int j = 0; j < mResX; ++j) {
				int index = i * mResX + j;
				Pixel& pixel = mPixelBuffer[index];
				pixel.color = image[index];
				pixel.weight = 1;
			}
		}
	}

	void SetFileName(const std::string &pFileName)
	{
		mFilename = pFileName;
	}

	void SaveImage();

private:

	void WriteToPixelBuffer() {
		for (int i = 0; i < mResY; ++i) {
			for (int j = 0; j < mResX; ++j) {
				int index = i * mResX + j;
				Pixel &pixel = mPixelBuffer[index];
				mImageBuffer[index] = pixel.color / pixel.weight + pixel.splat;
			}
		}
	}

};

GYT_NAMESPACE_END