#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <cassert>

struct STLTriangle {
    float normal[3];
    float v1[3];
    float v2[3];
    float v3[3];
};

bool equal(float a, float b, float epsilon = 1e-8) {

    return abs(a - b) < epsilon;
}

inline bool is_positive(float a) {
    return a > 0;
}

inline bool is_negative(float a) {
    return a < 0;
}

bool all_non_zero_same_sign(float a, float b, float c) {

    bool all_not_zero = (!equal(a, 0.)) && (!equal(b, 0.)) && (!equal(c, 0.));
    bool all_positive = is_positive(a) && is_positive(b) && is_positive(c);
    bool all_negative = is_negative(a) && is_negative(b) && is_negative(c);

    return all_not_zero && (all_positive || all_negative);
}

bool all_zero(float a, float b, float c) {
    return (equal(a, 0.)) && (equal(b, 0.)) && (equal(c, 0.));
}

bool non_zero_same_sign(float a, float b) {

    bool all_not_zero = (!equal(a, 0.)) && (!equal(b, 0.));
    bool all_positive = is_positive(a) && is_positive(b);
    bool all_negative = is_negative(a) && is_negative(b);

    return all_not_zero && (all_positive || all_negative);
}

float two_by_two_determinant(float a, float b, float c, float d) {
    /* Where the matrix = [[a, b],
    *                      [c, d]]
    */

   //std::cout << "----" << std::endl;
   //std::cout << a << ", " << b << std::endl;
   //std::cout << c << ", " << d << std::endl;
   //std::cout << "----" << std::endl;
   return (a * d) - (b * c);
}


struct CPUVector3 {

    public:
        CPUVector3(float, float, float);
        CPUVector3();
        float x, y, z;
        CPUVector3 operator-(const CPUVector3&) const;
        CPUVector3 operator-(CPUVector3&) const;
        CPUVector3 operator+(const CPUVector3&) const;
        CPUVector3 operator+(CPUVector3&) const;
        float operator*(CPUVector3&) const;
        float operator*(const CPUVector3&) const;
};


CPUVector3 cross_product(const CPUVector3& left, const CPUVector3& right) {

    float x = two_by_two_determinant(left.y, left.z, right.y, right.z);
    float y = -two_by_two_determinant(left.x, left.z, right.x, right.z);
    float z = two_by_two_determinant(left.x, left.y, right.x, right.y);

    return CPUVector3(x, y, z);
}

struct Triangle {

    public:
        Triangle();
        Triangle(const CPUVector3&, const CPUVector3&, const CPUVector3&);
        CPUVector3 vertex1;
        CPUVector3 vertex2;
        CPUVector3 vertex3;
};

Triangle::Triangle(): vertex1(CPUVector3(0, 0, 0)), vertex2(CPUVector3(0, 0, 0)), vertex3(CPUVector3(0, 0, 0)) {
}

Triangle::Triangle(const CPUVector3& vertex1, const CPUVector3& vertex2, const CPUVector3& vertex3): vertex1(vertex1), vertex2(vertex2), vertex3(vertex3) {
}

struct TriangleTriangleCollisionDetector {
    //https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tritri_tam.pdf

    public:
        TriangleTriangleCollisionDetector();
        TriangleTriangleCollisionDetector(const Triangle&, const Triangle&);
        bool check(void);
        const Triangle triangle_1;
        const Triangle triangle_2;
        CPUVector3 N1, N2;
        float d1, d2;
        float t1_d_v1, t1_d_v2, t1_d_v3;
        float t2_d_v1, t2_d_v2, t2_d_v3;
    
    private:
        CPUVector3 compute_signed_distances(const Triangle&, const CPUVector3&, float) const;
        CPUVector3 point_on_line(const CPUVector3&, float, const CPUVector3&, float) const;
        bool do_intervals_intersect(float, float, float, float);

};

CPUVector3 TriangleTriangleCollisionDetector::point_on_line(const CPUVector3& N1, float d1, const CPUVector3& N2, float d2) const {
    // Intersection of two planes: https://www.youtube.com/watch?v=O6O_64zIEYI
    // Symbolic matrix solver: https://www.symbolab.com

    float x, y, z;
    x = ((N1.y * d2) - (N2.y * d1)) / ((-N1.y * N2.x) + (N1.x * N2.y));
    y = ((-N1.x * d2) + (N2.x * d1)) / ((-N1.y * N2.x) + (N1.x * N2.y));
    z = 0.0;
    return CPUVector3(x, y, z);
}

bool TriangleTriangleCollisionDetector::do_intervals_intersect(float a, float b, float c, float d) {

    float lower1, upper1, lower2, upper2;
    lower1 = a <= b ? a : b;
    upper1 = b > a ? b : a;

    lower2 = c <= d ? c : d;
    upper2 = d > c ? d : c;

    if (lower1 >= lower2 && lower1 <= upper2) {
        return true;
    }
    else if (upper1 >= lower2 && upper1 <= upper2) {
        return true;
    }
    else if (lower2 >= lower1 && lower2 <= upper1) {
        return true;
    }
    else if (upper2 >= lower1 && upper2 <= upper1) {
        return true;
    }

    return false;
}

float line_segment_intersects_plane(CPUVector3 P, CPUVector3 Q, CPUVector3 N, float D) {

    // P = (2, 1, 3)
    /*
    P.x = 2.0;
    P.y = 1.0;
    P.z = 3.0;

    // Q = (5, 2, 1)
    Q.x = 5.0;
    Q.y = 2.0;
    Q.z = 1.0;

    // N = (1, -3, -5)
    N.x = 1.0;
    N.y = -3.0;
    N.z = -5.0;

    D = 4.0;
    */

   float numerator = -D + (N.x * P.x) + (N.y * P.y) + (N.z * P.z);
   float denominator = (N.x * (P.x - Q.x)) + (N.y * (P.y - Q.y)) + (N.z * (P.z - Q.z));
   return numerator / denominator;
}

bool TriangleTriangleCollisionDetector::check(void) {

    if (all_zero(t1_d_v1, t1_d_v2, t1_d_v3) || all_zero(t2_d_v1, t2_d_v2, t2_d_v3)) {
        // The triangles are co-planar. Need to add this check.
        //std::cout << "The triangles are coplanar." << std::endl;
        return false;
    }
    else if (all_non_zero_same_sign(t1_d_v1, t1_d_v2, t1_d_v3) || all_non_zero_same_sign(t2_d_v1, t2_d_v2, t2_d_v3)) {
        //std::cout << "Collision rules out by all the signed distances having the same signs." << std::endl;
        return false;
    }

    // Get two vertices from each triangle which lie on the same side.
    CPUVector3 T1_V1, T1_V2, T1_V3;
    if (non_zero_same_sign(t1_d_v1, t1_d_v2)) {
        T1_V1 = triangle_1.vertex1;
        T1_V2 = triangle_1.vertex2;
        T1_V3 = triangle_1.vertex3;
    }
    else if (non_zero_same_sign(t1_d_v1, t1_d_v3)) {
        T1_V1 = triangle_1.vertex1;
        T1_V2 = triangle_1.vertex3;
        T1_V3 = triangle_1.vertex2;
    }
    else {
        T1_V1 = triangle_1.vertex2;
        T1_V2 = triangle_1.vertex3;
        T1_V3 = triangle_1.vertex1;
    }

    CPUVector3 T2_V1, T2_V2, T2_V3;
    if (non_zero_same_sign(t2_d_v1, t2_d_v2)) {
        T2_V1 = triangle_2.vertex1;
        T2_V2 = triangle_2.vertex2;
        T2_V3 = triangle_2.vertex3;
    }
    else if (non_zero_same_sign(t2_d_v1, t2_d_v3)) {
        T2_V1 = triangle_2.vertex1;
        T2_V2 = triangle_2.vertex3;
        T2_V3 = triangle_2.vertex2;
    }
    else {
        T2_V1 = triangle_2.vertex2;
        T2_V2 = triangle_2.vertex3;
        T2_V3 = triangle_2.vertex1;
    }

    float t1_a = line_segment_intersects_plane(T1_V1, T1_V3, N2, d2);
    float t1_b = line_segment_intersects_plane(T1_V2, T1_V3, N2, d2);

    float t2_a = line_segment_intersects_plane(T2_V1, T2_V3, N1, d1);
    float t2_b = line_segment_intersects_plane(T2_V2, T2_V3, N1, d1);

    return do_intervals_intersect(t1_a, t1_b, t2_a, t2_b);
}

TriangleTriangleCollisionDetector::TriangleTriangleCollisionDetector(const Triangle& T1, const Triangle& T2): triangle_1(T1), triangle_2(T2) {

    N1 = cross_product((T1.vertex2 - T1.vertex1), (T1.vertex3 - T1.vertex1));
    N2 = cross_product((T2.vertex2 - T2.vertex1), (T2.vertex3 - T2.vertex1));

    d1 = -(N1 * T1.vertex1);
    d2 = -(N2 * T2.vertex1);

    CPUVector3 T0_vertex_distances = compute_signed_distances(T1, N2, d2);
    CPUVector3 T1_vertex_distances = compute_signed_distances(T2, N1, d1);

    t1_d_v1 = T0_vertex_distances.x;
    t1_d_v2 = T0_vertex_distances.y;
    t1_d_v3 = T0_vertex_distances.z;

    t2_d_v1 = T1_vertex_distances.x;
    t2_d_v2 = T1_vertex_distances.y;
    t2_d_v3 = T1_vertex_distances.z;
}

CPUVector3 TriangleTriangleCollisionDetector::compute_signed_distances(const Triangle& T, const CPUVector3& N, float d) const {

    float d1 = (N * T.vertex1) + d;
    float d2 = (N * T.vertex2) + d;
    float d3 = (N * T.vertex3) + d;

    return CPUVector3(d1, d2, d3);
}

CPUVector3::CPUVector3(float x, float y, float z): x(x), y(y), z(z) {
}

CPUVector3::CPUVector3(): x(0.), y(0.), z(0.) {
}

CPUVector3 CPUVector3::operator-(const CPUVector3& right) const {

    return CPUVector3(x - right.x, y - right.y, z - right.z);
}

CPUVector3 CPUVector3::operator+(const CPUVector3& right) const {

    return CPUVector3(x + right.x, y + right.y, z + right.z);
}

CPUVector3 CPUVector3::operator+(CPUVector3& right) const {

    return CPUVector3(x + right.x, y + right.y, z + right.z);
}

CPUVector3 CPUVector3::operator-(CPUVector3& right) const {

    return CPUVector3(x - right.x, y - right.y, z - right.z);
}

float CPUVector3::operator*(CPUVector3& right) const {

    return (x * right.x) + (y * right.y) + (z * right.z);
}

float CPUVector3::operator*(const CPUVector3& right) const {

    return (x * right.x) + (y * right.y) + (z * right.z);
}


class Mesh {

    public:
        Mesh();
        void load_from_file(std::string);
        void print_triangles(void);
        Triangle* triangles;
        //STLTriangle* stl_triangles;
        //Transformation
        int num_triangles;
        std::string file_name;
};

Mesh::Mesh(): triangles(nullptr), file_name(""), num_triangles(0) {
}

void Mesh::print_triangles(void) {

    if (!triangles) {
        std::cout << "There are no triangles in this mesh!" << std::endl;
        return;
    }

    std::cout << "Number of triangles: " << num_triangles << std::endl;
    for (int i = 0; i < num_triangles; i++) {
        auto& triangle = triangles[i];
        std::cout << "Vertex 1: " << triangle.vertex1.x << " " << triangle.vertex1.y << " " << triangle.vertex1.z << std::endl;
        std::cout << "Vertex 2: " << triangle.vertex2.x << " " << triangle.vertex2.y << " " << triangle.vertex2.z << std::endl;
        std::cout << "Vertex 3: " << triangle.vertex3.x << " " << triangle.vertex3.y << " " << triangle.vertex3.z << std::endl;
    }
}

std::vector<STLTriangle> parseBinarySTL(const std::string& filename) {

    std::vector<STLTriangle> triangles;
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return triangles;
    }

    // Read 80-byte header (ignore)
    file.seekg(80, std::ios::beg);

    // Read number of triangles (4-byte unsigned int)
    uint32_t numTriangles;
    file.read(reinterpret_cast<char*>(&numTriangles), sizeof(numTriangles));

    // Read each triangle
    for (uint32_t i = 0; i < numTriangles; ++i) {
        STLTriangle triangle;
        // Read normal vector (3 floats, 12 bytes)
        file.read(reinterpret_cast<char*>(&triangle.normal), sizeof(triangle.normal));
        // Read vertices (3 * 3 floats, 36 bytes)
        file.read(reinterpret_cast<char*>(&triangle.v1), sizeof(triangle.v1));
        file.read(reinterpret_cast<char*>(&triangle.v2), sizeof(triangle.v2));
        file.read(reinterpret_cast<char*>(&triangle.v3), sizeof(triangle.v3));
        // Read attribute count (2-byte unsigned int, ignored)
        file.seekg(2, std::ios::cur);
        triangles.push_back(triangle);
}

    return triangles;
}

void Mesh::load_from_file(std::string file_name) {

    std::vector<STLTriangle> read_triangles = parseBinarySTL(file_name);
    this->num_triangles = read_triangles.size();
    triangles = (Triangle*) malloc(sizeof(Triangle) * num_triangles);
    for (int i = 0; i < read_triangles.size(); i++) {
        auto v1 = CPUVector3(read_triangles[i].v1[0], read_triangles[i].v1[1], read_triangles[i].v1[2]);
        auto v2 = CPUVector3(read_triangles[i].v2[0], read_triangles[i].v2[1], read_triangles[i].v2[2]);
        auto v3 = CPUVector3(read_triangles[i].v3[0], read_triangles[i].v3[1], read_triangles[i].v3[2]);
        triangles[i] = Triangle(v1, v2, v3);
    }
}

void test_one() {

    Triangle T1 = Triangle(CPUVector3(1., 0.1, 0.), CPUVector3(0., 0.2, 0.), CPUVector3(0., 0.3, 1.));
    Triangle T2 = Triangle(CPUVector3(-10., -10., -10.1), CPUVector3(-200., -100., -100.1), CPUVector3(-330., -12.0, -1.5));
    TriangleTriangleCollisionDetector pair1 = TriangleTriangleCollisionDetector(T1, T2);

    Triangle T3 = Triangle(CPUVector3(1., 0.1, 0.), CPUVector3(-1., 0.2, 0.), CPUVector3(0., 0.3, 1.));
    Triangle T4 = Triangle(CPUVector3(0., 1., 0.1), CPUVector3(0., -1., 0.1), CPUVector3(0., 0.0, 1.5));
    TriangleTriangleCollisionDetector pair2 = TriangleTriangleCollisionDetector(T3, T4);

    assert(pair1.check() == 0);
    assert(pair2.check() != 0);
}

int check_meshes_for_collision(const Mesh& mesh1, const Mesh& mesh2) {

    int num_collisions = 0;
    for (int i = 0; i < mesh1.num_triangles; i++) {
        auto mesh_1_triangle = mesh1.triangles[i];

        for (int j = 0; j < mesh2.num_triangles; j++) {
            auto mesh_2_triangle = mesh2.triangles[j];
            mesh_2_triangle.vertex1.x += 10000;
            mesh_2_triangle.vertex1.y += 10000;
            mesh_2_triangle.vertex1.z += 10000;

            mesh_2_triangle.vertex2.x += 20000;
            mesh_2_triangle.vertex2.y += 20000;
            mesh_2_triangle.vertex2.z += 20000;

            mesh_2_triangle.vertex3.x += 30000;
            mesh_2_triangle.vertex3.y += 30000;
            mesh_2_triangle.vertex3.z += 30000;
            TriangleTriangleCollisionDetector detector = TriangleTriangleCollisionDetector(mesh_1_triangle, mesh_2_triangle);
            if (detector.check()) {
                num_collisions += 1;
            }
        }
    }

    return num_collisions;
}

void mesh_test() {

    Mesh fore_arm = Mesh();
    fore_arm.load_from_file("../meshes/forearm.stl");

    Mesh upper_arm = Mesh();
    upper_arm.load_from_file("../meshes/upperarm.stl");

    int num_collisions = check_meshes_for_collision(fore_arm, upper_arm);

    std::cout << "The number of collisions: " << num_collisions << std::endl;
}


int main() {

    mesh_test();

    return 0;
}