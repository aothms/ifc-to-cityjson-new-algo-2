#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <CGAL/Vector_3.h>
#include <CGAL/Polyhedron_3.h>

template <typename Kernel>
void createCube(CGAL::Polyhedron_3<Kernel>& P, double d) {
	typedef CGAL::Point_3<Kernel> Point;
	typedef std::array<CGAL::Point_3<Kernel>, 4> Quad;

	Quad bottom = { {
		Point(-d, -d, -d),
		Point(+d, -d, -d),
		Point(+d, +d, -d),
		Point(-d, +d, -d)
	} };

	Kernel::Vector_3 V(0, 0, d * 2);

	extrude(bottom, V, P);
}

template <typename Kernel>
void rotatePolyhedron(CGAL::Polyhedron_3<Kernel>& P, double angle, CGAL::Vector_3<Kernel> axis) {
	axis = axis / std::sqrt(CGAL::to_double(axis.squared_length()));

	double c = std::cos(angle);
	double s = std::sin(angle);

	// Compute the 4x4 rotation matrix
	CGAL::Aff_transformation_3<Kernel> rotation_matrix(
		c + (1 - c) * axis.x() * axis.x(),
		(1 - c) * axis.x() * axis.y() - s * axis.z(),
		(1 - c) * axis.x() * axis.z() + s * axis.y(),
		0.0,
		(1 - c) * axis.y() * axis.x() + s * axis.z(),
		c + (1 - c) * axis.y() * axis.y(),
		(1 - c) * axis.y() * axis.z() - s * axis.x(),
		0.0,
		(1 - c) * axis.z() * axis.x() - s * axis.y(),
		(1 - c) * axis.z() * axis.y() + s * axis.x(),
		c + (1 - c) * axis.z() * axis.z(),
		0.0
	);

	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = rotation_matrix.transform(it->point());
	}
}


template <typename Kernel>
void scalePolyhedron(CGAL::Polyhedron_3<Kernel>& P, typename Kernel::FT sx, typename Kernel::FT sy, typename Kernel::FT sz) {
	CGAL::Aff_transformation_3<Kernel> scale(
		sx, 0, 0, 0,
		0, sy, 0, 0,
		0, 0, sz, 0);
	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = scale.transform(it->point());
	}
}

template <typename Kernel>
void translatePolyhedron(CGAL::Polyhedron_3<Kernel>& P, CGAL::Vector_3<Kernel> v) {
	CGAL::Aff_transformation_3<Kernel> rotation_matrix(
		CGAL::TRANSLATION, v
	);

	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = rotation_matrix.transform(it->point());
	}
}


#endif