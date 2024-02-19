#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

struct STLTriangle {
    float normal[3];
    float v1[3];
    float v2[3];
    float v3[3];
    short attributeByteCount;
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

   return (a * d) - (b * c);
}


struct Vector3 {

    public:
        Vector3(float, float, float);
        Vector3();
        float x, y, z;
        Vector3 operator-(const Vector3&) const;
        Vector3 operator-(Vector3&) const;
        float operator*(Vector3&) const;
        float operator*(const Vector3&) const;
};


Vector3 cross_product(const Vector3& left, const Vector3& right) {
    
    float x = two_by_two_determinant(left.y, left.z, right.y, right.z);
    float y = two_by_two_determinant(left.x, left.z, right.x, right.z);
    float z = two_by_two_determinant(left.x, left.y, right.x, right.y);

    return Vector3(x, y, z);
}

struct Triangle {

    public:
        Triangle();
        Triangle(const Vector3&, const Vector3&, const Vector3&);
        Vector3 vertex1;
        Vector3 vertex2;
        Vector3 vertex3;
};

Triangle::Triangle(): vertex1(Vector3(0, 0, 0)), vertex2(Vector3(0, 0, 0)), vertex3(Vector3(0, 0, 0)) {
}

Triangle::Triangle(const Vector3& vertex1, const Vector3& vertex2, const Vector3& vertex3): vertex1(vertex1), vertex2(vertex2), vertex3(vertex3) {
}

struct TriangleTriangleCollisionDetector {
    //https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tritri_tam.pdf

    public:
        TriangleTriangleCollisionDetector();
        TriangleTriangleCollisionDetector(const Triangle&, const Triangle&);
        bool check(void);
        const Triangle triangle_1;
        const Triangle triangle_2;
        Vector3 N1, N2;
        float d1, d2;
        float t1_d_v1, t1_d_v2, t1_d_v3;
        float t2_d_v1, t2_d_v2, t2_d_v3;
    
    private:
        Vector3 compute_signed_distances(const Triangle&, const Vector3&, float) const;
        Vector3 point_on_line(const Vector3&, float, const Vector3&, float) const;
        bool do_intervals_intersect(float, float, float, float);

};

Vector3 TriangleTriangleCollisionDetector::point_on_line(const Vector3& N1, float d1, const Vector3& N2, float d2) const {
    // Intersection of two planes: https://www.youtube.com/watch?v=O6O_64zIEYI
    // Symbolic matrix solver: https://www.symbolab.com

    float x, y, z;
    x = ((N1.y * d2) - (N2.y * d1)) / ((-N1.y * N2.x) + (N1.x * N2.y));
    y = ((-N1.x * d2) + (N2.x * d1)) / ((-N1.y * N2.x) + (N1.x * N2.y));
    z = 0.0;
    return Vector3(x, y, z);
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

    return false;
}

bool TriangleTriangleCollisionDetector::check(void) {

    // Make sure the data has been populated.
    if (all_zero(t1_d_v1, t1_d_v2, t1_d_v3) || all_zero(t2_d_v1, t2_d_v2, t2_d_v3)) {
        // The triangles are co-planar. Need to add this check.
        std::cout << "The triangles are coplanar." << std::endl;
        return false;
    }
    else if (all_non_zero_same_sign(t1_d_v1, t1_d_v2, t1_d_v3) || all_non_zero_same_sign(t2_d_v1, t2_d_v2, t2_d_v3)) {
        std::cout << "Collision rules out by all the signed distances having the same signs." << std::endl;
        return false;
    }

    // Get two vertices from each triangle which lie on the same side.
    Vector3 T1_V1, T1_V2;
    float d_1_v1, d_1_v2;
    if (non_zero_same_sign(t1_d_v1, t1_d_v2)) {
        T1_V1 = triangle_1.vertex1;
        T1_V2 = triangle_1.vertex2;

        d_1_v1 = t1_d_v1;
        d_1_v2 = t1_d_v2;
    }
    else if (non_zero_same_sign(t1_d_v1, t1_d_v3)) {
        T1_V1 = triangle_1.vertex1;
        T1_V2 = triangle_1.vertex3;

        d_1_v1 = t1_d_v1;
        d_1_v2 = t1_d_v3;
    }
    else {
        T1_V1 = triangle_1.vertex2;
        T1_V2 = triangle_1.vertex3;

        d_1_v1 = t1_d_v2;
        d_1_v2 = t1_d_v3;
    }

    Vector3 T2_V1, T2_V2;
    float d_2_v1, d_2_v2;
    if (non_zero_same_sign(t2_d_v1, t2_d_v2)) {
        T2_V1 = triangle_2.vertex1;
        T2_V2 = triangle_2.vertex2;

        d_2_v1 = t2_d_v1;
        d_2_v2 = t2_d_v2;
    }
    else if (non_zero_same_sign(t2_d_v1, t2_d_v3)) {
        T2_V1 = triangle_2.vertex1;
        T2_V2 = triangle_2.vertex3;

        d_2_v1 = t2_d_v1;
        d_2_v2 = t2_d_v3;
    }
    else {
        T2_V1 = triangle_2.vertex2;
        T2_V2 = triangle_2.vertex3;

        d_2_v1 = t2_d_v2;
        d_2_v2 = t2_d_v3;
    }

    // Compute the intersection line.
    Vector3 D = cross_product(N1, N2); // The intersection line's direction.
    Vector3 O = point_on_line(N1, d1, N2, d2); // A point on the intersection line.
    
    float p_1_v1 = D * (T1_V1 - O);
    float p_1_v2 = D * (T1_V2 - O);

    float p_2_v1 = D * (T2_V1 - O);
    float p_2_v2 = D * (T2_V2 - O);

    float t1_a = p_1_v1 + ((p_1_v2 - p_1_v1) * (d_1_v1  / (d_1_v1 - d_1_v2)));
    float t1_b = p_1_v2 + ((p_1_v1 - p_1_v2) * (d_1_v2  / (d_1_v2 - d_1_v1)));

    float t2_a = p_2_v1 + ((p_2_v2 - p_2_v1) * (d_2_v1  / (d_2_v1 - d_2_v2)));
    float t2_b = p_2_v2 + ((p_2_v1 - p_2_v2) * (d_2_v2  / (d_2_v2 - d_2_v1)));

    return do_intervals_intersect(t1_a, t1_b, t2_a, t2_b);
}

TriangleTriangleCollisionDetector::TriangleTriangleCollisionDetector(const Triangle& T1, const Triangle& T2) {

    N1 = cross_product((T1.vertex2 - T1.vertex1), (T1.vertex3 - T1.vertex1));
    N2 = cross_product((T2.vertex2 - T2.vertex1), (T2.vertex3 - T2.vertex1));

    d1 = -(N1 * T1.vertex1);
    d2 = -(N2 * T2.vertex1);

    Vector3 T0_vertex_distances = compute_signed_distances(T1, N2, d2);
    Vector3 T1_vertex_distances = compute_signed_distances(T2, N1, d1);

    t1_d_v1 = T0_vertex_distances.x;
    t1_d_v2 = T0_vertex_distances.y;
    t1_d_v3 = T0_vertex_distances.z;

    t2_d_v1 = T1_vertex_distances.x;
    t2_d_v2 = T1_vertex_distances.y;
    t2_d_v3 = T1_vertex_distances.z;
}

Vector3 TriangleTriangleCollisionDetector::compute_signed_distances(const Triangle& T, const Vector3& N, float d) const {

    float d1 = (N * T.vertex1) + d;
    float d2 = (N * T.vertex2) + d;
    float d3 = (N * T.vertex3) + d;

    return Vector3(d1, d2, d3);
}

Vector3::Vector3(float x, float y, float z): x(x), y(y), z(z) {
}

Vector3::Vector3(): x(0.), y(0.), z(0.) {
}

Vector3 Vector3::operator-(const Vector3& right) const {

    return Vector3(x - right.x, y - right.y, z - right.z);
}

Vector3 Vector3::operator-(Vector3& right) const {

    return Vector3(x - right.x, y - right.y, z - right.z);
}

float Vector3::operator*(Vector3& right) const {

    return (x * right.x) + (y * right.y) + (z * right.z);
}

float Vector3::operator*(const Vector3& right) const {

    return (x * right.x) + (y * right.y) + (z * right.z);
}


class Mesh {

    public:
        Mesh();
        void load_from_file(std::string);
        void print_triangles(void);
        Triangle* triangles;
        STLTriangle* stl_triangles;
        //Transformation
        int num_triangles;
        std::string file_name;
};

Mesh::Mesh(): triangles(nullptr), file_name(""), num_triangles(0) {
}

void Mesh::print_triangles(void) {

    std::cout << "Number of triangles: " << num_triangles << std::endl;
    for (int i = 0; i < num_triangles; i++) {
        std::cout << "Normal: " << stl_triangles[i].normal[0] << " " << stl_triangles[i].normal[1] << " " << stl_triangles[i].normal[2] << std::endl;
        std::cout << "Vertex 1: " << stl_triangles[i].v1[0] << " " << stl_triangles[i].v1[1] << " " << stl_triangles[i].v1[2] << std::endl;
        std::cout << "Vertex 2: " << stl_triangles[i].v2[0] << " " << stl_triangles[i].v2[1] << " " << stl_triangles[i].v2[2] << std::endl;
        std::cout << "Vertex 3: " << stl_triangles[i].v3[0] << " " << stl_triangles[i].v3[1] << " " << stl_triangles[i].v3[2] << std::endl;
    }
}

void Mesh::load_from_file(std::string file_name) {

    std::ifstream file(file_name, std::ios::binary);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file.");
    }

    // Skip the header (80 bytes) in the binary STL file
    file.seekg(80, std::ios::beg);

    // Read the number of triangles
    uint32_t num_triangles;
    file.read(reinterpret_cast<char*>(&num_triangles), sizeof(uint32_t));

    // Allocate space for the triangles.
    stl_triangles = (STLTriangle*) malloc(sizeof(STLTriangle) * num_triangles);

    // Read each triangle
    for (int i = 0; i < num_triangles; ++i) {
        file.read(reinterpret_cast<char*>(&stl_triangles[i]), sizeof(STLTriangle));
    }

    this->num_triangles = num_triangles;

    // Close the file
    file.close();
}

int main() {
    
    Mesh mesh = Mesh();
    mesh.load_from_file("./forearm.stl");
    // mesh.print_triangles();

    // Vertices of triangle 1.
    auto t1_v1 = mesh.stl_triangles[0].v1;
    auto t1_v2 = mesh.stl_triangles[0].v2;
    auto t1_v3 = mesh.stl_triangles[0].v3;

    // Vertices of triangle 2.
    auto t2_v1 = mesh.stl_triangles[1].v1;
    auto t2_v2 = mesh.stl_triangles[1].v2;
    auto t2_v3 = mesh.stl_triangles[1].v3;

    std::cout << "[" << t1_v1[0] << ", " << t1_v1[1] << ", " << t1_v1[2] << "]" << std::endl;

    Triangle T1 = Triangle(Vector3(12., 24., 36.), Vector3(45., 54., 62.), Vector3(71., 89., 99.));
    Triangle T2 = Triangle(Vector3(-13., -24., -35.), Vector3(-46., -57., -68.), Vector3(-68., -69., -70.));
    
    Triangle T3 = Triangle(Vector3(1., 0., 0.), Vector3(-1., 0., 0.), Vector3(0., 0., 1.));
    Triangle T4 = Triangle(Vector3(0., 3., 0.), Vector3(0., -3., 0.), Vector3(0., 0., 0.5));


    TriangleTriangleCollisionDetector info1 = TriangleTriangleCollisionDetector(T1, T2);
    std::cout << info1.d1 << std::endl;
    auto result = info1.check();
    std::cout << "Collision: " << result << std::endl;

    TriangleTriangleCollisionDetector info2 = TriangleTriangleCollisionDetector(T3, T4);
    std::cout << info2.d1 << std::endl;
    result = info2.check();
    std::cout << "Collision: " << result << std::endl;

    return 0;
}