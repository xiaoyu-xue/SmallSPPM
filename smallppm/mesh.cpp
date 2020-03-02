#include "mesh.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

NAMESPACE_BEGIN

void Mesh::LoadFromFile(std::string inputfile, bool reverseVertices) {
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());

	if (!warn.empty()) {
		std::cout << warn << std::endl;
	}
	if (!err.empty()) {
		std::cerr << err << std::endl;
	}
	if (!ret) {
		exit(1);
	}
	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			int fv = shapes[s].mesh.num_face_vertices[f];

			// Loop over vertices in the face.
			int i = (int)vertices.size(), j = i + 1, k = i + 2;
			if (reverseVertices) {
				std::swap(j, k);
			}
			bool has_normal = false;
			for (int v = 0; v < fv; v++) {
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				float vx = attrib.vertices[3 * idx.vertex_index + 0];
				float vy = attrib.vertices[3 * idx.vertex_index + 1];
				float vz = attrib.vertices[3 * idx.vertex_index + 2];
				vertices.push_back(Vec3(vx, vy, vz));
				if (idx.normal_index != -1) {
					float nx = attrib.normals[3 * idx.normal_index + 0];
					float ny = attrib.normals[3 * idx.normal_index + 1];
					float nz = attrib.normals[3 * idx.normal_index + 2];
					has_normal = true;
					normals.push_back(Vec3(nx, ny, nz));
				}
				else {
					// should be consistent
					assert(!has_normal);
					has_normal = false;
				}
				float tx = 0.0f, ty = 0.0f;
				if (idx.texcoord_index != -1) {
					tx = attrib.texcoords[2 * idx.texcoord_index + 0];
					ty = attrib.texcoords[2 * idx.texcoord_index + 1];
				}
				uvs.push_back(Vector2(tx, ty));
			}
			if (!has_normal) {
				Vec3* a = &vertices[vertices.size() - 1] - 2;
				Vec3 generated_normal = Cross(a[1] - a[0], a[2] - a[0]);
				if (generated_normal.Length() > 1e-6f) {
					generated_normal = generated_normal.Norm();
				}
				for (int v = 0; v < fv; v++) {
					normals.push_back(generated_normal);
				}
			}

			untransformedTriangles.push_back(
				Triangle(nullptr, nullptr,
					vertices[i], vertices[j], vertices[k],
					normals[i], normals[j], normals[k],
					uvs[i], uvs[j], uvs[k]));
			index_offset += fv;
		}
	}
}

void Mesh::SetMaterial(const std::shared_ptr<Material> material) {
	this->material = material;
}
void Mesh::SetEmission(const Vec3 emission) {
	this->isEmission = true;
	this->emission = emission;
}

NAMESPACE_END