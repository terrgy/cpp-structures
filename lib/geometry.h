#include <vector>
#include <cmath>
#include <set>
#include <typeinfo>


static bool doubleComparison(double a, double b) {
    static const double EPS = 1e-9;
    return fabs(a - b) < EPS;
}

static std::pair<double, double> angleToCosSin(double angle) {
    static const int PI_TO_DEGREE = 180;
    angle = angle * M_PI / PI_TO_DEGREE;
    return {cos(angle), sin(angle)};
}

struct Vector;
class Line;

struct Point {
    double x, y;

    Point(double x, double y) : x(x), y(y) {}
    Point() : Point(0, 0) {}

    bool operator==(const Point& another) const {
        return doubleComparison(x, another.x) && doubleComparison(y, another.y);
    }

    bool operator<(const Point& another) const {
        static const double EPS = 1e-9;
        return (x - another.x < -EPS) || (doubleComparison(x, another.x) && (y - another.y < -EPS));
    }

    void rotate(const Point& center, double angle);

    void rotate(const Point& center, double cos_, double sin_);

    void translateByVector(const Vector& vector);

    void reflect(const Line& axis);

    void scale(const Point& center, double coefficient);

    void reflect(const Point& center) {
        scale(center, -1);
    }

    static double sqrDistance(const Point& point1, const Point& point2) {
        return (point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y);
    }

    static double distance(const Point& point1, const Point& point2) {
        return sqrt(sqrDistance(point1, point2));
    }

    static Point getMiddlePoint(const Point& point1, const Point& point2) {
        static const double HALF_SCALE = 0.5;
        Point middle = point2;
        middle.scale(point1, HALF_SCALE);
        return middle;
    }
};

struct Vector : public Point {
public:
    explicit Vector(double x, double y) : Point(x, y) {}
    explicit Vector(const Point& point) : Vector(point.x, point.y) {}
    Vector(const Point& start, const Point& end) : Vector(end.x - start.x, end.y - start.y) {}

    double sqrLength() const {
        return x * x + y * y;
    }

    double length() const {
        return sqrt(sqrLength());
    }

    void setLength(double new_length) {
        double ratio = new_length / length();
        *this *= ratio;

    }

    void rotate(double angle) {
        auto [cos_, sin_] = angleToCosSin(angle);
        rotate(cos_, sin_);
    }

    void rotate(double cos_, double sin_) {
        double tmp_x = x;
        x = x * cos_ - y * sin_;
        y = tmp_x * sin_ + y * cos_;
    }

    Vector& operator+=(const Vector& another) {
        x += another.x;
        y += another.y;
        return *this;
    }

    Vector& operator*=(double coefficient) {
        x *= coefficient;
        y *= coefficient;
        return *this;
    }

    Vector& operator/=(double coefficient) {
        x /= coefficient;
        y /= coefficient;
        return *this;
    }

    Vector operator*(double coefficient) {
        Vector result(*this);
        result *= coefficient;
        return result;
    }

    Vector operator/(double coefficient) {
        Vector result(*this);
        result /= coefficient;
        return result;
    }

    static double crossProduct(const Vector& vector1, const Vector& vector2) {
        return vector1.x * vector2.y - vector1.y * vector2.x;
    }

    static double scalarProduct(const Vector& vector1, const Vector& vector2) {
        return vector1.x * vector2.x + vector1.y * vector2.y;
    }

    static bool isVectorBetweenVectors(const Vector& target, const Vector& lower, const Vector& upper) {
        double lower_cross_product = crossProduct(lower, target);
        return ((lower_cross_product > 0) && (crossProduct(upper, target) < 0)) ||
               (doubleComparison(lower_cross_product, 0) && (scalarProduct(lower, target) > 0));
    }

    static bool areAnglesBetweenVectorsEqual(const std::pair<Vector, Vector>& angle1,
                                             const std::pair<Vector, Vector>& angle2) {
        return doubleComparison(
                crossProduct(angle1.first, angle1.second) * angle2.first.length() * angle2.second.length(),
                crossProduct(angle2.first, angle2.second) * angle1.first.length() * angle1.second.length());
    }

    static Vector getMiddleVector(const Vector& vector1, const Vector& vector2) {
        Vector middle = vector1;
        middle.setLength(vector2.length());
        middle += vector2;
        return middle;
    }
};

void Point::translateByVector(const Vector &vector) {
    x += vector.x;
    y += vector.y;
}

void Point::rotate(const Point &center, double angle) {
    auto [cos_, sin_] = angleToCosSin(angle);
    Point::rotate(center, cos_, sin_);
}

void Point::rotate(const Point& center, double cos_, double sin_) {
    Vector point_vector(center, *this);
    point_vector.rotate(cos_, sin_);
    *this = center;
    translateByVector(point_vector);
}

void Point::scale(const Point &center, double coefficient) {
    Vector point_vector(center, *this);
    point_vector *= coefficient;
    *this = center;
    translateByVector(point_vector);
}

class Line {
public:
    double a, b, c;

    Line(double a_coeff, double b_coeff, double c_coeff) : a(a_coeff), b(b_coeff), c(c_coeff) {}
    Line(const Point& point1, const Point& point2) :
            a(point1.y - point2.y), b(point2.x - point1.x), c(-a * point1.x - b * point1.y) {}
    Line(double ang, double shift) : a(-ang), b(1), c(-shift) {}
    Line(const Point& point, double ang) : Line(ang, point.y - ang * point.x) {}

    bool operator==(const Line& another) const {
        return doubleComparison(a * another.b, another.a * b) && doubleComparison(c * another.b, another.c * b) &&
               doubleComparison(c * another.a, another.c * a);
    }

    bool containsPoint(const Point& point) const {
        return doubleComparison(a * point.x + b * point.y + c, 0);
    }

    Vector getPerpendicularVector() const {
        return Vector(a, b);
    }

    Vector getParallelVector() const {
        return Vector(b, -a);
    }

    static double distance(const Line& line, const Point& point) {
        return fabs(line.a * point.x + line.b * point.y + line.c) / sqrt(line.a * line.a + line.b * line.b);
    }

    static Point intersection(const Line& line1, const Line& line2)
    {
        double x = (line2.c * line1.b - line1.c * line2.b) / (line1.a * line2.b - line2.a * line1.b);
        double y = (line1.a * line2.c - line1.c * line2.a) / (line1.b * line2.a - line1.a * line2.b);
        return {x, y};
    }

    static Line getPerpendicularBisector(const Point& point1, const Point& point2) {
        return getLineByPerpendicularVector(Point::getMiddlePoint(point1, point2), Vector(point1, point2));
    }

    static Line getBisector(const Point& from, const Point& left, const Point& right) {
        Point end = from;
        end.translateByVector(Vector::getMiddleVector(
                Vector(from, left), Vector(from, right)));
        return {from, end};
    }

    static Line getMedian(const Point& from, const Point& to1, const Point& to2) {
        return {from, Point::getMiddlePoint(to1, to2)};
    }

    static Point getAltitudePoint(const Point& from, const Line& line) {
        Line altitude = getLineByPerpendicularVector(from, line.getParallelVector());
        return intersection(line, altitude);
    }

    static Line getAltitude(const Point& from, const Line& line) {
        return {from, getAltitudePoint(from, line)};
    }

    static Line getLineByPerpendicularVector(const Point& point, const Vector& perpendicular_vector) {
        return {perpendicular_vector.x, perpendicular_vector.y,
                -perpendicular_vector.x * point.x - perpendicular_vector.y * point.y};
    }

    static bool isPointOnSegment(const Point& point, const Point& left, const Point& right) {
        if (!Line(left, right).containsPoint(point)) {
            return false;
        }
        return Vector::scalarProduct(Vector(left, point), Vector(right, point)) < 0;
    }
};

void Point::reflect(const Line& axis) {
    Point altitude_point = Line::getAltitudePoint(*this, axis);
    altitude_point.scale(*this, 2);
    *this = altitude_point;
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool equals(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;

    virtual ~Shape() = default;
};

bool operator==(const Shape& shape1, const Shape& shape2) {
    return shape1.equals(shape2);
}

class Polygon : public Shape {
private:
    std::vector<Point> vertices;

    static size_t _nextIndexStatic(size_t current, size_t size) {
        return current + 1 < size ? current + 1 : 0;
    }

    size_t _nextIndex(size_t current) const {
        return _nextIndexStatic(current, verticesCount());
    }

    static size_t _prevIndexStatic(size_t current, size_t size) {
        return current ? current - 1 : size - 1;
    }

    size_t _prevIndex(size_t current) const {
        return _prevIndexStatic(current, verticesCount());
    }

    std::pair<Vector, Vector> _getAngleVectors(size_t prev, size_t current, size_t next) const {
        return {Vector(vertices[current], vertices[prev]), Vector(vertices[current], vertices[next])};
    }

    static void _moveindices(size_t& prev, size_t& current, size_t (*nextIndexFunc)(size_t, size_t), size_t size) {
        prev = nextIndexFunc(prev, size);
        current = nextIndexFunc(current, size);
    }

    static void _moveindices(size_t& prev, size_t& current, size_t& next, size_t (*nextIndexFunc)(size_t, size_t), size_t size) {
        _moveindices(prev, current, nextIndexFunc, size);
        next = nextIndexFunc(next, size);
    }

    bool _checkAnglesSimilarity(const Polygon& another, size_t shift, size_t (*nextIndexFunc)(size_t, size_t)) const {
        size_t prev = shift, current = nextIndexFunc(prev, verticesCount()), next = nextIndexFunc(current, verticesCount());
        size_t an_prev = 0, an_current = _nextIndex(an_prev), an_next = _nextIndex(an_current);
        for (size_t i = 0; i < verticesCount(); ++i) {
            std::pair<Vector, Vector> angle1 = _getAngleVectors(prev, current, next);
            std::pair<Vector, Vector> angle2 = another._getAngleVectors(an_prev, an_current, an_next);
            if (!Vector::areAnglesBetweenVectorsEqual(angle1, angle2)) {
                return false;
            }
            _moveindices(prev, current, next, nextIndexFunc, verticesCount());
            _moveindices(an_prev, an_current, an_next, _nextIndexStatic, another.verticesCount());
        }
        return true;
    }

    bool _checkSidesSimilarity(const Polygon& another, size_t shift, size_t (*nextIndexFunc)(size_t, size_t)) const {
        size_t prev = shift, current = nextIndexFunc(prev, verticesCount());
        size_t an_prev = 0, an_current = _nextIndex(an_prev);
        for (size_t i = 0; i < verticesCount(); ++i) {
            double sqr_dist1 = Point::sqrDistance(vertices[prev], vertices[current]);
            double sqr_dist2 = Point::sqrDistance(another.vertices[an_prev], another.vertices[an_current]);
            if (!doubleComparison(sqr_dist1, sqr_dist2)) {
                return false;
            }
            _moveindices(prev, current,  nextIndexFunc, verticesCount());
            _moveindices(an_prev, an_current,  _nextIndexStatic, another.verticesCount());
        }
        return true;
    }

    bool _equalDirection(const Polygon& another, size_t shift, size_t (*nextIndexFunc)(size_t, size_t)) const {
        for (size_t i = 0, current = shift; i < verticesCount(); ++i, current = nextIndexFunc(current, verticesCount())) {
            if (vertices[current] != another.vertices[i]) {
                return false;
            }
        }
        return true;
    }

    bool _isCongruentOneSideCheck(const Polygon& another, size_t (*nextIndexFunc)(size_t, size_t)) const {
        for (size_t shift = 0; shift < verticesCount(); ++shift) {
            if (_checkAnglesSimilarity(another, shift, nextIndexFunc) &&
                _checkSidesSimilarity(another, shift, nextIndexFunc))
            {
                return true;
            }
        }
        return false;
    }

    bool _isCongruentOnePolygonCheck(const Polygon& another) const {
        return _isCongruentOneSideCheck(another, _nextIndexStatic) ||
               _isCongruentOneSideCheck(another, _prevIndexStatic);
    }

public:
    Polygon() = default;
    explicit Polygon(const std::vector<Point>& vertices) : vertices(vertices) {}
    Polygon(std::initializer_list<Point> lst) : vertices(lst) {}
    template <typename... Points>
    explicit Polygon(Points... points) : vertices({points...}) {}

    size_t verticesCount() const {
        return vertices.size();
    }

    const std::vector<Point>& getVertices() const {
        return vertices;
    }

    double perimeter() const final {
        double perimeter = 0;
        for (size_t i = 0; i < verticesCount(); ++i) {
            perimeter += Point::distance(vertices[i],
                                         vertices[_nextIndex(i)]);
        }
        return perimeter;
    }

    double area() const final {
        double area = 0;
        for (size_t i = 0; i < verticesCount(); ++i)
        {
            area += Vector::crossProduct(Vector(vertices[i]),
                                         Vector(vertices[_nextIndex(i)]));
        }
        return fabs(area) / 2;
    }


    bool equals(const Polygon& another) const {
        for (size_t shift = 0; shift < verticesCount(); ++shift) {
            if (_equalDirection(another, shift, _nextIndexStatic) || _equalDirection(another, shift, _prevIndexStatic)) {
                return true;
            }
        }
        return false;
    }

    bool equals(const Shape& another) const final {
        const auto* polygon = dynamic_cast<const Polygon*>(&another);
        if (!polygon) {
            return false;
        }
        return equals(*polygon);
    }

    bool isConvex() const {
        bool rotate_left = false, rotate_right = false;
        for (size_t current = 0, next = 1, next_next = 2;
             current < verticesCount();
             ++current, next = _nextIndex(next), next_next = _nextIndex(next_next))
        {
            if (Vector::crossProduct(
                    Vector(vertices[current], vertices[next]),
                    Vector(vertices[next], vertices[next_next])) < 0)
            {
                rotate_left = true;
            } else {
                rotate_right = true;
            }
        }
        return rotate_left xor rotate_right;
    }

    bool isCongruentTo(const Polygon& another) const {
        if (verticesCount() != another.verticesCount()) {
            return false;
        }
        Polygon reversed = *this;
        reversed.reflect(Line(Point(0, 0), Point(1, 0)));
        return _isCongruentOnePolygonCheck(another) || reversed._isCongruentOnePolygonCheck(another);
    }

    bool isSimilarTo(const Polygon& another) const {
        if (verticesCount() != another.verticesCount()) {
            return false;
        }
        Polygon another_copy = another;
        another_copy.scale(Point(0, 0), perimeter() / another.perimeter());
        return isCongruentTo(another_copy);
    }

    bool isCongruentTo(const Shape& another) const final {
        try {
            return isCongruentTo(dynamic_cast<const Polygon&>(another));
        } catch (std::bad_cast&) {
            return false;
        }
    }

    bool isSimilarTo(const Shape& another) const final {
        try {
            return isSimilarTo(dynamic_cast<const Polygon&>(another));
        } catch (std::bad_cast&) {
            return false;
        }
    }

    bool containsPoint(const Point& point) const final {
        static const double EPS = 1e-5;
        for (auto vertex : vertices) {
            if (vertex == point) {
                return true;
            }
        }

        for (size_t i = 0; i < verticesCount(); ++i) {
            if (Line::isPointOnSegment(point, vertices[i], vertices[_nextIndex(i)])) {
                return true;
            }
        }

        bool intersections_flag = false;
        Polygon rotated_copy = *this;
        rotated_copy.rotate(point, EPS);
        for (size_t i = 0, j = verticesCount() - 1; i < verticesCount(); j = i++) {
            if (((rotated_copy.vertices[i].y >= point.y - EPS * EPS) !=
                 (rotated_copy.vertices[j].y >= point.y - EPS * EPS)) &&
                (point.x - EPS <= (rotated_copy.vertices[j].x - rotated_copy.vertices[i].x) *
                                  (point.y - rotated_copy.vertices[i].y) /
                                  (rotated_copy.vertices[j].y - rotated_copy.vertices[i].y) +
                                  rotated_copy.vertices[i].x))
            {
                intersections_flag = !intersections_flag;
            }
        }
        return intersections_flag;
    }

    void rotate(const Point& center, double angle) final {
        auto [cos_, sin_] = angleToCosSin(angle);
        for (size_t i = 0; i < verticesCount(); ++i) {
            vertices[i].rotate(center, cos_, sin_);
        }
    }

    void reflect(const Point& center) final {
        for (size_t i = 0; i < verticesCount(); ++i) {
            vertices[i].reflect(center);
        }
    }

    void reflect(const Line& axis) final {
        for (size_t i = 0; i < verticesCount(); ++i) {
            vertices[i].reflect(axis);
        }
    }

    void scale(const Point& center, double coefficient) final {
        for (size_t i = 0; i < verticesCount(); ++i) {
            vertices[i].scale(center, coefficient);
        }
    }
};

class Ellipse : public Shape {
private:
    Point focus1, focus2;
    double a, c;

public:
    Ellipse(const Point& focus1, const Point& focus2, double dist_to_focuses) :
            focus1(focus1), focus2(focus2), a(dist_to_focuses / 2), c(Point::distance(focus1, focus2) / 2) {}

    std::pair<Point, Point> focuses() const {
        return {focus1, focus2};
    }

    std::pair<Line, Line> directrices() const {
        Vector translate_vector(focus1, focus2);
        translate_vector.setLength(a / eccentricity() - c);
        Point directrix_point1(focus1), directrix_point2(focus2);
        directrix_point1.translateByVector(translate_vector * -1);
        directrix_point2.translateByVector(translate_vector);
        return {Line::getLineByPerpendicularVector(directrix_point1, translate_vector),
                Line::getLineByPerpendicularVector(directrix_point2, translate_vector)};
    }

    double eccentricity() const {
        return c / a;
    }

    Point center() const {
        return Point::getMiddlePoint(focus1, focus2);
    }

    double getSemiMajorAxis() const {
        return a;
    }

    double getSemiMinorAxis() const {
        return sqrt(a * a - c * c);
    }

    double perimeter() const final {
        static const int COEFFICIENT = 3;
        double semi_minor_axis = getSemiMinorAxis();
        return M_PI * (COEFFICIENT * (a + semi_minor_axis) -
                       sqrt((COEFFICIENT * a + semi_minor_axis) * (a + COEFFICIENT * semi_minor_axis)));
    }

    double area() const final {
        return M_PI * a * getSemiMinorAxis();
    }

    bool equals(const Ellipse& another) const {
        return (std::min(focus1, focus2) == std::min(another.focus1, another.focus2)) &&
               (std::max(focus1, focus2) == std::max(another.focus1, another.focus2)) && doubleComparison(a, another.a);
    }

    bool equals(const Shape& another) const final {
        const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
        if (!ellipse) {
            return false;
        }
        return equals(*ellipse);
    }

    bool isCongruentTo(const Ellipse& another) const {
        return doubleComparison(a, another.a) && doubleComparison(c, another.c);
    }

    bool isCongruentTo(const Shape& another) const final {
        try{
            return isCongruentTo(dynamic_cast<const Ellipse&>(another));
        } catch (std::bad_cast&) {
            return false;
        }
    }

    bool isSimilarTo(const Ellipse& another) const {
        return doubleComparison(eccentricity(), another.eccentricity());
    }

    bool isSimilarTo(const Shape& another) const final {
        try{
            return isSimilarTo(dynamic_cast<const Ellipse&>(another));
        } catch (std::bad_cast&) {
            return false;
        }
    }

    bool containsPoint(const Point& point) const final {
        return Point::distance(focus1, point) + Point::distance(focus2, point) <= 2 * a;
    }

    void rotate(const Point& center, double angle) final {
        auto [cos_, sin_] = angleToCosSin(angle);
        focus1.rotate(center, cos_, sin_);
        focus2.rotate(center, cos_, sin_);
    }

    void reflect(const Point& center) final {
        focus1.reflect(center);
        focus2.reflect(center);
    }

    void reflect(const Line& axis) final {
        focus1.reflect(axis);
        focus2.reflect(axis);
    }

    void scale(const Point& center, double coefficient) final {
        focus1.scale(center, coefficient);
        focus2.scale(center, coefficient);
        a *= coefficient;
    }
};

class Circle : public Ellipse {
public:
    Circle(const Point& center, double radius) : Ellipse(center, center, 2 * radius) {}

    double radius() const {
        return getSemiMajorAxis();
    }
};

class Rectangle : public Polygon {
private:
    static Point getRectanglePointLeftDiagonal(const Point& point1, const Point& point2, double ratio) {
        double diag_sqr_length = Point::sqrDistance(point1, point2);
        double b = diag_sqr_length / (ratio * ratio + 1), a = diag_sqr_length - b;
        if (a > b) {
            std::swap(a, b);
        }
        double kx = 2 * (point2.x - point1.x);
        double ky = 2 * (point2.y - point1.y);
        double w = a - b - point1.x * point1.x + point2.x * point2.x - point1.y * point1.y + point2.y * point2.y;

        double q = 1 + kx * kx / (ky * ky);
        double r = 2 * (point1.y * kx / ky - w * kx / (ky * ky) - point1.x);
        double s = point1.x * point1.x + point1.y * point1.y - a + w * w / (ky * ky) - 2 * point1.y * w / ky;

        static const int DISCRIMINANT_COEFFICIENT = 4;
        double D = sqrt(r * r - DISCRIMINANT_COEFFICIENT * q * s);

        double x1 = (-r + D) / (2 * q), x2 = (-r - D) / (2 * q);

        double y1 = w / ky - kx * x1 / ky, y2 = w / ky - kx * x2 / ky;

        Point candidate1(x1, y1), candidate2(x2, y2);
        return Vector::crossProduct(Vector(point1, point2), Vector(point1, candidate1)) >= 0 ? candidate1 : candidate2;
    }

public:
    Rectangle(const Point& point1, const Point& point2, double ratio) :
            Polygon{point1, getRectanglePointLeftDiagonal(point1, point2, ratio),
                    point2, getRectanglePointLeftDiagonal(point2, point1, ratio)} {}

    Point center() const {
        Vector point_translate(getVertices()[0], getVertices()[2]);
        point_translate /= 2;
        Point center = getVertices()[0];
        center.translateByVector(point_translate);
        return center;
    }

    std::pair<Line, Line> diagonals() const {
        static const std::pair<size_t, size_t> first_indices = {0, 2};
        static const std::pair<size_t, size_t> second_indices = {1, 3};
        return {Line(getVertices()[first_indices.first], getVertices()[first_indices.second]),
                Line(getVertices()[second_indices.first], getVertices()[second_indices.second])};
    }
};

class Square : public Rectangle {
public:
    Square(const Point& point1, const Point& point2) : Rectangle(point1, point2, 1) {}

    double getSideLength() const {
        return Point::distance(getVertices()[0], getVertices()[1]);
    }

    double getDiagonalLength() const {
        return Point::distance(getVertices()[0], getVertices()[2]);
    }

    Circle circumscribedCircle() const {
        return {center(), getSideLength() / 2};
    }

    Circle inscribedCircle() const {
        return {center(), getDiagonalLength() / 2};
    }
};

class Triangle : public Polygon {
public:
    Triangle(const Point& point1, const Point& point2, const Point& point3) : Polygon(point1, point2, point3) {}

    Circle circumscribedCircle() const {
        Line line1 = Line::getPerpendicularBisector(getVertices()[0], getVertices()[1]);
        Line line2 = Line::getPerpendicularBisector(getVertices()[1], getVertices()[2]);
        Point circle_center = Line::intersection(line1, line2);
        return {circle_center, Point::distance(circle_center, getVertices()[0])};
    }

    Circle inscribedCircle() const {
        Line line1 = Line::getBisector(getVertices()[0], getVertices()[1], getVertices()[2]);
        Line line2 = Line::getBisector(getVertices()[1], getVertices()[0], getVertices()[2]);
        Point circle_center = Line::intersection(line1, line2);
        Point altitude = Line::getAltitudePoint(circle_center, Line(getVertices()[0], getVertices()[1]));
        return {circle_center, Point::distance(circle_center, altitude)};
    }

    Point centroid() const {
        Line line1 = Line::getMedian(getVertices()[0], getVertices()[1], getVertices()[2]);
        Line line2 = Line::getMedian(getVertices()[1], getVertices()[0], getVertices()[2]);
        return Line::intersection(line1, line2);
    }

    Point orthocenter() const {
        Line line1 = Line::getAltitude(getVertices()[0], Line(getVertices()[1], getVertices()[2]));
        Line line2 = Line::getAltitude(getVertices()[1], Line(getVertices()[0], getVertices()[2]));
        return Line::intersection(line1, line2);
    }

    Line EulerLine() const {
        return {circumscribedCircle().center(), orthocenter()};
    }

    Circle ninePointsCircle() const {
        Point circle_center = Point::getMiddlePoint(circumscribedCircle().center(), orthocenter());
        return {circle_center, Point::distance(circle_center,
                                               Point::getMiddlePoint(getVertices()[0], getVertices()[1]))};
    }
};