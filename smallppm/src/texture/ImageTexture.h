#pragma once
#include "Texture.h"
#include "math/MathUtils.h"
#include "image/ImageIO.h"
#include "math/Linagl.h"

NAMESPACE_BEGIN

enum class WrapMode {
	REPEAT, BLACK, CLAMP
};

enum class FilterMode {
	NEAREST,
	BILINEAR
};

template<typename T>
class ImageTexture : public Texture<T> {
public:
	ImageTexture(const Array2D<T>& image, WrapMode wrapMode = WrapMode::CLAMP, FilterMode filterMode = FilterMode::BILINEAR) :
		texture(image), wrapMode(wrapMode), filterMode(filterMode) {
	}

	ImageTexture(const std::string& filename, WrapMode wrapMode = WrapMode::CLAMP, FilterMode filterMode = FilterMode::BILINEAR) :
		wrapMode(wrapMode), filterMode(filterMode) {
		texture = ImageIO::LoadTexture(filename);
	}

	T Sample(const Vec2& coord) const override {
		int resX = texture.res[1];
		int resY = texture.res[0];
		if (filterMode == FilterMode::BILINEAR) {
			return Triangle(coord[0], coord[1]);
		}
		else {
			return Nearest(coord[0], coord[1]);
		}
	}

	T Sample(const Intersection& isect) const override {
		return Sample(isect.uv);
	}

	T Texel(int u, int v) const {
		if (wrapMode == WrapMode::CLAMP) {
			u = Clamp(u, 0, texture.res[1] - 1);
			v = Clamp(v, 0, texture.res[0] - 1);
		}
		return texture[v][u];
	}

	T Triangle(float u, float v) const {
		int resX = texture.res[1];
		int resY = texture.res[0];
		float s = resX * u - 0.5f;
		float t = resY * v - 0.5f;
		int s0 = std::floor(s), t0 = std::floor(t);
		float ds = s - s0, dt = t - t0;
		T texel = Lerp(dt,
			Lerp(ds, Texel(s0, t0), Texel(s0 + 1, t0)),
			Lerp(ds, Texel(s0, t0 + 1), Texel(s0 + 1, t0 + 1)));
		return texel;
	}

	T Nearest(float u, float v) const {
		int resX = texture.res[1];
		int resY = texture.res[0];
		T texel = Texel(int(u * resX), int(v * resY));
		return texel;
	}

	int Height() const {
		return texture.res[0];
	}

	int Width() const {
		return texture.res[1];
	}

	Vec3 ElementAt(int i, int j) const {
		return texture[texture.res[0] - j - 1][i];
		//return texture[j][i];

	}
public:
	Array2D<T> texture;
	const WrapMode wrapMode;
	const FilterMode filterMode;
};

NAMESPACE_END