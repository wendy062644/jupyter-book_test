#include <bits/stdc++.h>
using namespace std;

int file_num = 1;

struct Point {
    double x, y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y || fabs(x - other.x) < 0.01 && fabs(y - other.y) < 0.01;
    }

    Point& operator=(const Point& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    bool operator<(const Point& other) const {
        return x < other.x || x == other.x && y < other.y;
    }

    bool operator!=(const Point& other) const {
        return fabs(x - other.x) >= 0.01 && fabs(y - other.y) >= 0.01;
    }
};

struct Edge {
    Point start, end, A, B;
};

vector<Point> points, convexhull_point;
vector<Edge> edges;
int Max_num = 100000, Min_num = -100000;

pair<Point, Point> calculatePerpendicularBisector(const Point& A, const Point& B) {
    double x = (A.x + B.x) / 2;
    double y = (A.y + B.y) / 2;

    Point startPoint, endPoint;

    const double Min_num = -100000;
    const double Max_num = 100000;

    // 水平
    if (B.y - A.y == 0) {
        startPoint.x = x;
        startPoint.y = Min_num;
        endPoint.x = x;
        endPoint.y = Max_num;
        return {startPoint, endPoint};
    }

    // 垂直
    if (B.x - A.x == 0) {
        startPoint.x = Min_num;
        startPoint.y = y;
        endPoint.x = Max_num;
        endPoint.y = y;
        return {startPoint, endPoint};
    }

    double m = -(B.x - A.x) / (B.y - A.y);

    if (m > 0) {
        if ((Max_num - y) / m + x <= Max_num) {
            endPoint.x = x + (Max_num - y) / m;
            endPoint.y = Max_num;
        } else {
            endPoint.x = Max_num;
            endPoint.y = y + m * (Max_num - x);
        }
        if ((y - Min_num) / m + x >= Min_num) {
            startPoint.x = x - (y - Min_num) / m;
            startPoint.y = Min_num;
        } else {
            startPoint.x = Min_num;
            startPoint.y = y - m * (x - Min_num);
        }
    } else {
        if ((y - Min_num) / m + x <= Max_num) {
            endPoint.x = x - (y - Min_num) / m;
            endPoint.y = Min_num;
        } else {
            endPoint.x = Max_num;
            endPoint.y = y + m * (Max_num - x);
        }

        if ((Max_num - y) / m + x >= Min_num) {
            startPoint.x = x + (Max_num - y) / m;
            startPoint.y = Max_num;
        } else {
            startPoint.x = Min_num;
            startPoint.y = y - m * (x - Min_num);
        }
    }

    return {startPoint, endPoint};
}

Point calculateMidPoint(const Point& A, const Point& B) {
    return {(A.x + B.x) / 2.0, (A.y + B.y) / 2.0};
}

Point calculateIntersection(double slope1, Point mid1, double slope2, Point mid2) {
    if (isinf(slope1)) {
        return {mid1.x, slope2 * (mid1.x - mid2.x) + mid2.y};
    }
    if (isinf(slope2)) {
        return {mid2.x, slope1 * (mid2.x - mid1.x) + mid1.y};
    }

    // y - mid1.y = slope1 * (x - mid1.x)
    // y - mid2.y = slope2 * (x - mid2.x)
    double x = (slope2 * mid2.x - slope1 * mid1.x + mid1.y - mid2.y) / (slope2 - slope1);
    double y = slope1 * (x - mid1.x) + mid1.y;
    return {x, y};
}

pair<Point, Point> extendLine(const Point& a, const Point& b) {
    double m;
    if (b.x != a.x) {
        m = (b.y - a.y) / (b.x - a.x);
    } else {
        m = numeric_limits<double>::infinity();
    }

    Point newA = a;
    Point newB = b;

    if (b.x != a.x) {
        newA.x = -10000;
        newA.y = a.y + m * (newA.x - a.x);
    } else { // 垂直
        newA.y = -10000;
    }

    if (b.x != a.x) {
        newB.x = 10000;
        newB.y = a.y + m * (newB.x - a.x);
    } else { // 垂直
        newB.y = 10000;
    }

    return {newA, newB};
}

double crossProduct(const Point& p, const Point& a, const Point& b) {
    // 水平
    if (a.y == b.y) {
        if (p.y < a.y) return 1;
        if (p.y > a.y) return -1;
        return 0;
    }

    // 垂直
    if (a.x == b.x) {
        if (p.x < a.x) return 1;
        if (p.x > a.x) return -1;
        return 0;
    }

    double computedX = (a.x - b.x) / (a.y - b.y) * (p.y - b.y) + b.x;

    if (p.x < computedX) {
        return 1;
    } else if (p.x > computedX) {
        return -1;
    } else {
        return 0;
    }
}


// 判斷點 P 是否在線段 AB 上
bool isOnSegment(const Point& A, const Point& B, const Point& P) {

    if(fabs(A.x - B.x) < 0.0001 && abs(P.x - A.x) < 0.0001) return true;
    if(fabs(A.y - B.y) < 0.0001 && abs(P.y - A.y) < 0.0001) return true;

    double crossProduct = (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);

    if (fabs(crossProduct) > 0.0001) {
        return false;
    }

    return min(A.x, B.x) <= P.x && P.x <= max(A.x, B.x) &&
           min(A.y, B.y) <= P.y && P.y <= max(A.y, B.y);
}

pair<Point, Point> adjustPointToBoundary(const Point& p1, const Point& p2) {
    Point startPoint = p1, endPoint = p2;
    if (p2.y - p1.y == 0) { // 水平
        startPoint.x = 0;
        startPoint.y = p2.y;
        endPoint.x = Max_num;
        endPoint.y = p2.y;
        return {startPoint, endPoint};
    }
    if (p2.x - p1.x == 0) { // 垂直
        startPoint.x = p2.x;
        startPoint.y = 0;
        endPoint.x = p2.x;
        endPoint.y = Max_num;
        return {startPoint, endPoint};
    }
    double m = (p1.y-p2.y) / (p1.x-p2.x);

    if (p2.x > Max_num) {  // 右邊界
        endPoint.x = Max_num;
        endPoint.y = p1.y + (Max_num - p1.x) * m;
    } else if (p2.x < 0) {  // 左邊界
        endPoint.x = 0;
        endPoint.y = p1.y + (0 - p1.x) * m;
    }

    if (endPoint.y > Max_num) {  // 上邊界
        endPoint.y = Max_num;
        endPoint.x = p1.x + (Max_num - p1.y) / m;
    } else if (endPoint.y < 0) {  // 下邊界
        endPoint.y = 0;
        endPoint.x = p1.x - (p1.y / m);
    }
    if (p1.x > Max_num) {
        startPoint.x = Max_num;
        startPoint.y = p2.y + (Max_num - p2.x) * m;
    } else if (p1.x < 0) {
        startPoint.x = 0;
        startPoint.y = p2.y + (0 - p2.x) * m;
    }

    if (startPoint.y > Max_num) {
        startPoint.y = Max_num;
        startPoint.x = p2.x + (Max_num - p2.y) / m;
    } else if (startPoint.y < 0) {
        startPoint.y = 0;
        startPoint.x = p2.x - (p2.y / m);
    }
    return {startPoint, endPoint};
}

vector<Point> convexHull(Point* begin, Point* end) {
    vector<Point> points(begin, end + 1);

    if (points.size() < 4) return points;

    sort(points.begin(), points.end());

    vector<Point> hull;

    for (const auto& point : points) {
        while (hull.size() >= 2 && crossProduct(hull[hull.size() - 2], hull[hull.size() - 1], point) <= 0) {
            hull.pop_back();
        }
        hull.push_back(point);
    }

    size_t t = hull.size();
    for (int i = points.size() - 2; i >= 0; --i) {
        while (hull.size() > t && crossProduct(hull[hull.size() - 2], hull[hull.size() - 1], points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }

    hull.pop_back();

    return hull;
}

double distanceToSegment(const Point& A, const Point& B, const Point& P) {
    double ABx = B.x - A.x;
    double ABy = B.y - A.y;
    double APx = P.x - A.x;
    double APy = P.y - A.y;

    double dotProduct = ABx * APx + ABy * APy;
    double lengthSquared = ABx * ABx + ABy * ABy;

    double t = max(0.0, min(1.0, dotProduct / lengthSquared));
    double projX = A.x + t * ABx;
    double projY = A.y + t * ABy;

    return sqrt((P.x - projX) * (P.x - projX) + (P.y - projY) * (P.y - projY));
}

void Swap(Point &a, Point &b) {
    if(a.y > b.y) swap(a, b);
}

vector<Edge> two_point(const Point A, const Point B) {
    pair<Point, Point> t = calculatePerpendicularBisector(A, B);
    Swap(t.first, t.second);
    Point a = A, b = B;
    Swap(a, b);
    return vector<Edge>{Edge{t.first, t.second, a, b}};
}

vector<Edge> trangle(const Point A, const Point B, const Point C, Point &t) {
    if(B.x == A.x && C.x == A.x
        || B.y == A.y && C.y == A.y
        || (B.y - A.y)*(C.x - B.x) == (C.y - B.y)*(B.x - A.x)) { // 垂直 or 水平
        return vector<Edge>{
            Edge{two_point(A, B)[0].start, two_point(A, B)[0].end, A, B},
            Edge{two_point(B, C)[0].start, two_point(B, C)[0].end, B, C}
        };
    }
    double slope1 = -(B.x - A.x) / (B.y - A.y);
    double slope2 = -(C.x - B.x) / (C.y - B.y);

    Point mid1 = calculateMidPoint(A, B);
    Point mid2 = calculateMidPoint(B, C);

    pair<Point, Point> bisector1 = calculatePerpendicularBisector(A, B);

    t = calculateIntersection(slope1, mid1, slope2, mid2);
    if(t.x >= 0 && t.x <= 600 && t.y >= 0 && t.y <= 600) {
        if(distanceToSegment(bisector1.first, C, t) > distanceToSegment(bisector1.second, C, t)) {
            bisector1.first = t;
        }else{
            bisector1.second = t;
        }
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, A, B}
            //Edge{bisector2.first, bisector2.second, b, c},
            //Edge{bisector3.first, bisector3.second, a, c}
        };
    }
    return vector<Edge>{};
}

Point IntersectionPoint(const Edge& A, const Edge& B) {
    // 判斷是否為垂直線或水平線
    bool A_vertical = fabs(A.end.x - A.start.x) < 0.01;
    bool A_horizontal = fabs(A.end.y - A.start.y) < 0.01;
    bool B_vertical = fabs(B.end.x - B.start.x) < 0.01;
    bool B_horizontal = fabs(B.end.y - B.start.y) < 0.01;

    double x, y;

    if (A_vertical && B_vertical) {
        // 都垂直
        return Point{-12345, -12345};
        throw invalid_argument("Both edges are vertical and parallel, no intersection.");
    } else if (A_horizontal && B_horizontal) {
        // 都水平
        return Point{-12345, -12345};
        throw invalid_argument("Both edges are horizontal and parallel, no intersection.");
    } else if (A_vertical) {
        x = A.start.x;
        if (B_horizontal) {
            y = B.start.y;
        } else {
            double slope_B = (B.end.y - B.start.y) / (B.end.x - B.start.x);
            y = slope_B * (x - B.start.x) + B.start.y;
        }
    } else if (B_vertical) {
        x = B.start.x;
        if (A_horizontal) {
            y = A.start.y;
        } else {
            double slope1 = (A.end.y - A.start.y) / (A.end.x - A.start.x);
            y = slope1 * (x - A.start.x) + A.start.y;
        }
    } else if (A_horizontal) {
        y = A.start.y;
        double slope2 = (B.end.y - B.start.y) / (B.end.x - B.start.x);
        x = (y - B.start.y) / slope2 + B.start.x;
    } else if (B_horizontal) {
        y = B.start.y;
        double slope1 = (A.end.y - A.start.y) / (A.end.x - A.start.x);
        x = (y - A.start.y) / slope1 + A.start.x;
    } else {
        double slope1 = (A.end.y - A.start.y) / (A.end.x - A.start.x);
        double slope2 = (B.end.y - B.start.y) / (B.end.x - B.start.x);

        if (slope1 == slope2) {
            return {-12345, -12345};
            throw invalid_argument("Edges are parallel, no intersection.");
        }

        x = ((slope1 * A.start.x - A.start.y) - (slope2 * B.start.x - B.start.y)) / (slope1 - slope2);
        y = slope1 * (x - A.start.x) + A.start.y;
    }

    return Point{x, y};
}

pair<Point, Point> HyperPlane(const vector<Point>& rightpoint, const vector<Point>& leftpoint) {
    pair<Point, Point> best_line;
    double max_y = -numeric_limits<double>::infinity();
    double best_slope = numeric_limits<double>::infinity();

    for (const auto& lp : leftpoint) {
        for (const auto& rp : rightpoint) {
            double m = (rp.y - lp.y) / (rp.x - lp.x);
            double b = lp.y - m * lp.x;

            bool is_valid = true;
            for (const auto& pt : leftpoint) {
                if (pt.y > m * pt.x + b) {
                    is_valid = false;
                    break;
                }
            }
            for (const auto& pt : rightpoint) {
                if (pt.y > m * pt.x + b) {
                    is_valid = false;
                    break;
                }
            }

            if (is_valid) {
                double line_y = max(lp.y, rp.y);

                if (line_y > max_y ||
                    (line_y == max_y && m == best_slope &&
                     (rp.x < best_line.first.x || (rp.x == best_line.first.x && lp.x > best_line.second.x)))) {
                    max_y = line_y;
                    best_slope = m;
                    best_line = {rp, lp};
                }
            }
        }
    }
    return best_line;
}

optional<Point> doLinesIntersect(const Edge& A, const Edge& B) {
    // 計算每個點與線段的相對位置
    double d1 = crossProduct(B.start, B.end, A.start);
    double d2 = crossProduct(B.start, B.end, A.end);
    double d3 = crossProduct(A.start, A.end, B.start);
    double d4 = crossProduct(A.start, A.end, B.end);

    // 判斷是否相交
    if (d1 * d2 < 0 && d3 * d4 < 0) {
        // 計算交點
        double a1 = A.end.y - A.start.y;
        double b1 = A.start.x - A.end.x;
        double c1 = a1 * A.start.x + b1 * A.start.y;

        double a2 = B.end.y - B.start.y;
        double b2 = B.start.x - B.end.x;
        double c2 = a2 * B.start.x + b2 * B.start.y;

        double determinant = a1 * b2 - a2 * b1;

        if (determinant != 0) { // 兩條線段相交
            double x = (b2 * c1 - b1 * c2) / determinant;
            double y = (a1 * c2 - a2 * c1) / determinant;
            Point intersection = {x, y};

            // 確保交點在兩條線段的範圍內
            if (isOnSegment(A.start, A.end, intersection) && isOnSegment(B.start, B.end, intersection)) {
                return intersection;
            }
        }
    }

    // 判斷是否在端點上
    if (d1 == 0 && isOnSegment(B.start, B.end, A.start)) return A.start;
    if (d2 == 0 && isOnSegment(B.start, B.end, A.end)) return A.end;
    if (d3 == 0 && isOnSegment(A.start, A.end, B.start)) return B.start;
    if (d4 == 0 && isOnSegment(A.start, A.end, B.end)) return B.end;

    return nullopt;
}

vector<Edge> constructConvexHull(vector<Edge>& voronoiEdges, bool Left) {
    vector<Edge> Convexhull;

    for(int i = 0; i < voronoiEdges.size(); i++) {
        if(voronoiEdges[i].start.x == voronoiEdges[i].end.x || voronoiEdges[i].start.y == voronoiEdges[i].end.y) {
            Convexhull.push_back(voronoiEdges[i]);
            //voronoiEdges.erase(voronoiEdges.begin() + i--);
            continue;
        }
        if(voronoiEdges[i].start.y > voronoiEdges[i].end.y) {
            swap(voronoiEdges[i].start, voronoiEdges[i].end);
        }
        double m = (voronoiEdges[i].end.y - voronoiEdges[i].start.y)/(voronoiEdges[i].end.x - voronoiEdges[i].start.x);
        if(voronoiEdges[i].start.y == 0 || voronoiEdges[i].start.x == 0 || voronoiEdges[i].start.y == 600 || voronoiEdges[i].start.x == 600 ||
           voronoiEdges[i].end.y == 0 || voronoiEdges[i].end.x == 0 || voronoiEdges[i].end.y == 600 || voronoiEdges[i].end.x == 600) {
            Convexhull.push_back(voronoiEdges[i]);
            //voronoiEdges.erase(voronoiEdges.begin() + i--);
        }
    }
    return Convexhull;
}

vector<Point> ConvexHull_Point(const vector<Edge>& voronoiEdges) {
    auto compare = [](const Point& p1, const Point& p2) {
        return p1.y > p2.y || p1.y == p2.y && p1.x > p2.x;
    };
    set<Point, decltype(compare)> ConvexHull(compare);
    for(Edge i : voronoiEdges) {
        ConvexHull.insert(i.A);
        ConvexHull.insert(i.B);
    }
    vector<Point> result;
    for(Point i : ConvexHull) {
        result.push_back(i);
    }
    return result;
}

double calculateslope(Edge E) {
    return (E.start.y - E.end.y)/(E.start.x - E.end.x);
}

vector<Edge> recursiveVoronoi(int L, int R) {
    convexhull_point.clear();
    if(L+1 == R) return vector<Edge>{};
    vector<Edge> voronoi;
    if (R - L == 2) {
        pair<Point, Point> t = calculatePerpendicularBisector(points[L], points[L+1]);
        Swap(t.first, t.second);
        Point A = points[L], B = points[L+1];
        Swap(A, B);
        return vector<Edge>{Edge{t.first, t.second, A, B}};
    } else if (R - L == 3) {
        if(points[L+1].x == points[L].x && points[L+2].x == points[L].x
           || points[L+1].y == points[L].y && points[L+2].y == points[L].y
           || (points[L+1].y - points[L].y)*(points[L+2].x - points[L+1].x) == (points[L+2].y - points[L+1].y)*(points[L+1].x - points[L].x)) { // 垂直 or 水平
            return vector<Edge>{
                Edge{recursiveVoronoi(L, L+2)[0].start, recursiveVoronoi(L, L+2)[0].end, points[L], points[L+1]},
                Edge{recursiveVoronoi(L+1, R)[0].start, recursiveVoronoi(L+1, R)[0].end, points[L+1], points[L+2]}
            };
        }
        double slope1 = -(points[L+1].x - points[L].x) / (points[L+1].y - points[L].y);
        double slope2 = -(points[L+2].x - points[L+1].x) / (points[L+2].y - points[L+1].y);

        Point mid1 = calculateMidPoint(points[L], points[L+1]);
        Point mid2 = calculateMidPoint(points[L+1], points[L+2]);

        pair<Point, Point> bisector1 = calculatePerpendicularBisector(points[L], points[L+1]);
        pair<Point, Point> bisector2 = calculatePerpendicularBisector(points[L+1], points[L+2]);
        pair<Point, Point> bisector3 = calculatePerpendicularBisector(points[L], points[L+2]);

        Point intersection = calculateIntersection(slope1, mid1, slope2, mid2);

        if(intersection.x >= 0 && intersection.x <= 600 && intersection.y >= 0 && intersection.y <= 600) {
            double AB2 = pow(points[L].x - points[L+1].x, 2) + pow(points[L].y - points[L+1].y, 2);
            double BC2 = pow(points[L+1].x - points[L+2].x, 2) + pow(points[L+1].y - points[L+2].y, 2);
            double CA2 = pow(points[L].x - points[L+2].x, 2) + pow(points[L].y - points[L+2].y, 2);
            if(AB2 >= BC2 && AB2 >= CA2) {
                swap(bisector1, bisector3);
                swap(points[L+1], points[L+2]);
            }
            else if(BC2 >= AB2 && BC2 >= CA2) {
                swap(bisector2, bisector3);
                swap(points[L], points[L+1]);
            }
            Point t = calculateMidPoint(points[L], points[L+1]);
            if(isOnSegment(intersection, bisector1.first, t)) {
                bisector1.second = intersection;
            }
            else {
                bisector1.first = intersection;
            }

            t = calculateMidPoint(points[L+1], points[L+2]);
            if(isOnSegment(intersection, bisector2.first, t)) {
                bisector2.second = intersection;
            }
            else {
                bisector2.first = intersection;
            }

            Swap(bisector3.second, bisector3.first);
            if(AB2 + BC2 == CA2 || AB2 + CA2 == BC2 || BC2 + CA2 == AB2) { // 直角
                if(calculateslope(Edge{points[L+1], intersection}) > 0) {
                    if(points[L+1].x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }
                else {
                    if(points[L+1].x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
            }
            else if(AB2 + BC2 < CA2 || AB2 + CA2 < BC2 || BC2 + CA2 < AB2) { // 鈍角
                if(calculateslope(Edge{calculateMidPoint(points[L], points[L+2]), intersection}) > 0) {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }
                else {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
            }
            else { // 銳角
                if(calculateslope(Edge{calculateMidPoint(points[L], points[L+2]), intersection}) > 0) {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
                else {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }

            }
        }
        else {
            sort(points.begin()+L, points.begin()+L+3, [](Point a, Point b) {
                    return a.y > b.y;
                 });
            return vector<Edge>{
                Edge{recursiveVoronoi(L, L+2)[0].start, recursiveVoronoi(L, L+2)[0].end, points[L+1], points[L]},
                Edge{recursiveVoronoi(L+1, R)[0].start, recursiveVoronoi(L+1, R)[0].end, points[L+2], points[L+1]}
            };
        }

        Swap(bisector1.first, bisector1.second);
        Swap(bisector2.first, bisector2.second);
        Swap(bisector3.first, bisector3.second);
        Point A = points[L], B = points[L+1], C = points[L+2];
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, (A.y < B.y ? A:B), (A.y > B.y ? A:B)},
            Edge{bisector2.first, bisector2.second, (B.y < C.y ? B:C), (B.y > C.y ? B:C)},
            Edge{bisector3.first, bisector3.second, (A.y < C.y ? A:C), (A.y > C.y ? A:C)}
        };
    }
    int mid = (L + R) / 2;
    vector<Edge> RightEdge = recursiveVoronoi(L, mid);
    vector<Edge> LeftEdge = recursiveVoronoi(mid, R);

    vector<Edge> RightConvexhull_Edge = constructConvexHull(RightEdge, 0);
    vector<Edge> LeftConvexhull_Edge = constructConvexHull(LeftEdge, 1);

    vector<Point> RightConvexhull = ConvexHull_Point(RightEdge);
    vector<Point> LeftConvexhull = ConvexHull_Point(LeftEdge);

    convexhull_point.insert(convexhull_point.end(), RightConvexhull.begin(), RightConvexhull.end());
    convexhull_point.insert(convexhull_point.end(), LeftConvexhull.begin(), LeftConvexhull.end());

    vector<Point> R_Points, L_Points;

    ofstream outfile("convexhull/convexhull" + to_string(file_num) + ".txt");
    for (const auto& edge : RightConvexhull) {
        outfile << "R " << edge.x << " " << edge.y << "\n";
    }
    outfile << "\n";
    for (const auto& edge : LeftConvexhull) {
        outfile << "L " << edge.x << " " << edge.y << "\n";
    }
    outfile << "\n";
    outfile.close();

    vector<Edge> mid_edge, temp;

    pair<Point, Point> hyperplane = HyperPlane(LeftConvexhull, RightConvexhull);
    Point A = hyperplane.first, B = hyperplane.second, P;

    L_Points.push_back(hyperplane.first);
    R_Points.push_back(hyperplane.second);
    LeftConvexhull.erase(find(LeftConvexhull.begin(), LeftConvexhull.end(), hyperplane.first));
    RightConvexhull.erase(find(RightConvexhull.begin(), RightConvexhull.end(), hyperplane.second));

    while(LeftConvexhull.size() || RightConvexhull.size()) {
        pair<Point, Point> hyp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);
        Swap(hyp.first, hyp.second);
        if(mid_edge.size()) {
            hyp.second = prev(mid_edge.end()) -> start;
        }
        Point r = {0, -1}, l = {-1, -1};
        int r_index = -1, l_index = -1, no_outline = 0;
        for(int i = 0; i < RightEdge.size(); i++) {
            Point intersect = IntersectionPoint(Edge{hyp.first, hyp.second}, RightEdge[i]);
            if(intersect.y <= hyp.second.y && isOnSegment(RightEdge[i].start, RightEdge[i].end, intersect) && intersect.y > r.y) {
                r = intersect;
                r_index = i;
            }
            else if(intersect.y < hyp.second.y && intersect.x >= 0 && intersect.y >= 0 && intersect.x <= 600 && intersect.y <= 600
                && intersect.y >= r.y && isOnSegment(RightEdge[i].start, RightEdge[i].end, intersect)) {
                r = intersect;
                r_index = i;
            }
        }
        for(int i = 0; i < LeftEdge.size(); i++) {
            Point intersect = IntersectionPoint(Edge{hyp.first, hyp.second}, LeftEdge[i]);
            if(intersect.y <= hyp.second.y && intersect.y > l.y && isOnSegment(LeftEdge[i].start, LeftEdge[i].end, intersect)) {
                l = intersect;
                l_index = i;
            }
            else if(intersect.y < hyp.second.y && intersect.x >= 0 && intersect.y >= 0 && intersect.x <= 600 && intersect.y <= 600
                && intersect.y >= l.y && isOnSegment(LeftEdge[i].start, LeftEdge[i].end, intersect)) {
                l = intersect;
                l_index = i;
            }
        }
        if(r.y == -1 && l.y == -1) break;
        if(r.y > l.y) {
            if(r.y < 600) {
                if(crossProduct(RightEdge[r_index].end, hyp.first, hyp.second) < 0) {
                    P = RightEdge[r_index].end;
                    RightEdge[r_index].end = r;
                }
                else {
                    P = RightEdge[r_index].start;
                    RightEdge[r_index].start = r;
                }
                hyp.first = r;
                no_outline = 1;
            }
            temp.push_back(RightEdge[r_index]);
            RightEdge.erase(RightEdge.begin() + r_index);
            for(int i = 0; i < RightConvexhull.size(); i++) {
                pair<Point, Point> t = calculatePerpendicularBisector(hyperplane.first, RightConvexhull[i]);
                Point p = IntersectionPoint(Edge{hyp.first, hyp.second}, Edge{t.first, t.second});
                if(p == r) {
                    hyperplane.second = RightConvexhull[i];
                    R_Points.push_back(RightConvexhull[i]);
                    RightConvexhull.erase(RightConvexhull.begin() + i);
                    break;
                }
            }
        }
        else if(l.y > r.y) {
            if(l.y < 600) {
                if(crossProduct(LeftEdge[l_index].end, hyp.first, hyp.second) < 0) {
                    P = LeftEdge[l_index].start;
                    LeftEdge[l_index].start = l;
                }
                else {
                    P = LeftEdge[l_index].end;
                    LeftEdge[l_index].end = l;
                }
                hyp.first = l;
                no_outline = 1;
            }
            temp.push_back(LeftEdge[l_index]);
            LeftEdge.erase(LeftEdge.begin() + l_index);
            for(int i = 0; i < LeftConvexhull.size(); i++) {
                pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull[i], hyperplane.second);
                Point p = IntersectionPoint(Edge{hyp.first, hyp.second}, Edge{t.first, t.second});
                if(p == l) {
                    hyperplane.first = LeftConvexhull[i];
                    L_Points.push_back(LeftConvexhull[i]);
                    LeftConvexhull.erase(LeftConvexhull.begin() + i--);
                    break;
                }
            }
        }
        else {
            if(crossProduct(RightEdge[r_index].end, hyp.first, hyp.second) > 0)
                RightEdge[r_index].end = r;
            else
                RightEdge[r_index].start = r;

            if(crossProduct(LeftEdge[l_index].end, hyp.first, hyp.second) < 0)
                LeftEdge[l_index].end = l;
            else
                LeftEdge[l_index].start = l;

            temp.push_back(RightEdge[r_index]);
            temp.push_back(LeftEdge[l_index]);
            RightEdge.erase(RightEdge.begin() + r_index);
            LeftEdge.erase(LeftEdge.begin() + l_index);
            bool b = 0;
            for(int i = 0; i < RightConvexhull.size() && !b; i++) {
                for(int j = 0; j < LeftConvexhull.size(); j++) {
                    pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull[j], RightConvexhull[i]);
                    Point p = IntersectionPoint(Edge{hyp.first, hyp.second}, Edge{t.first, t.second});
                    if(p == r || p.x == -12345 && p.y == -12345) {
                        hyperplane.first = LeftConvexhull[j];
                        L_Points.push_back(LeftConvexhull[i]);
                        LeftConvexhull.erase(LeftConvexhull.begin() + j);
                        hyperplane.second = RightConvexhull[i];
                        R_Points.push_back(RightConvexhull[i]);
                        RightConvexhull.erase(RightConvexhull.begin() + i);
                        b = 1;
                        break;
                    }
                }

            }
        }
        for(int i = 0; i < LeftConvexhull.size(); i++) {
            if(LeftConvexhull[i].y > hyperplane.first.y) {
                L_Points.push_back(LeftConvexhull[i]);
                LeftConvexhull.erase(LeftConvexhull.begin() + i--);
            }
        }
        for(int i = 0; i < RightConvexhull.size(); i++) {
            if(RightConvexhull[i].y > hyperplane.second.y) {
                R_Points.push_back(RightConvexhull[i]);
                RightConvexhull.erase(RightConvexhull.begin() + i--);
            }
        }
        if(no_outline) {
            mid_edge.push_back(Edge{hyp.first, hyp.second, A, B});
        }
    }

    pair<Point, Point> hyp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);
    Swap(hyp.first, hyp.second);
    if(mid_edge.size()) {
        hyp.second = prev(mid_edge.end()) -> start;
    }

    for(int i = 0; i < temp.size(); i++) {
        Point intersect = IntersectionPoint(Edge{hyp.first, hyp.second}, temp[i]);
        if(isOnSegment(hyp.first, hyp.second, intersect) && isOnSegment(temp[i].start, temp[i].end, intersect)
           && fabs(intersect.x - hyp.second.x) > 0.01 && intersect.y <= hyp.second.y) {
            hyp.first = intersect;
            if(intersect.x >= 0 && intersect.x <= 600 && intersect.y >= 0 && intersect.y <= 600) {
                if(crossProduct(temp[i].B, hyp.first, hyp.second) > 0) {
                    if(crossProduct(temp[i].end, hyp.first, hyp.second) > 0) {
                        temp[i].start = intersect;
                    }
                    else {
                        temp[i].end = intersect;
                    }
                }
                else {
                    if(crossProduct(temp[i].end, hyp.first, hyp.second) < 0) {
                        temp[i].start = intersect;
                    }
                    else {
                        temp[i].end = intersect;
                    }
                }
            }
        }
    }
    mid_edge.push_back(Edge{hyp.first, hyp.second, hyperplane.first, hyperplane.second});

    for(int i = 0; i < R_Points.size(); i++) {
        hyp = calculatePerpendicularBisector(hyperplane.first, R_Points[i]);
        Swap(hyp.first, hyp.second);
        if(mid_edge.size()) {
            hyp.second = prev(mid_edge.end()) -> start;
        }
        Point intersect = IntersectionPoint(Edge{hyp.first, hyp.second}, *prev(mid_edge.end()));
        if(intersect == prev(mid_edge.end()) -> start) {
            hyp.second = intersect;
            mid_edge.push_back(Edge{hyp.first, hyp.second, hyperplane.first, R_Points[i]});
            break;
        }
    }

    for(int i = 0; i < L_Points.size(); i++) {
        hyp = calculatePerpendicularBisector(L_Points[i], hyperplane.second);
        Swap(hyp.first, hyp.second);
        if(mid_edge.size()) {
            hyp.second = prev(mid_edge.end()) -> start;
        }
        Point intersect = IntersectionPoint(Edge{hyp.first, hyp.second}, *prev(mid_edge.end()));
        if(intersect == prev(mid_edge.end()) -> start) {
            hyp.second = intersect;
            mid_edge.push_back(Edge{hyp.first, hyp.second, L_Points[i], hyperplane.second});
            break;
        }
    }

    outfile.open("hyperplane/hyperplane" + to_string(file_num++) + ".txt");
    for (const auto& edge : mid_edge) {
        outfile << edge.end.x << " " << edge.end.y << "\n";
    }

    outfile << mid_edge[mid_edge.size()-1].start.x << " " << mid_edge[mid_edge.size()-1].start.y << "\n";
    outfile.close();

    voronoi.insert(voronoi.end(), RightEdge.begin(), RightEdge.end());
    voronoi.insert(voronoi.end(), LeftEdge.begin(), LeftEdge.end());
    voronoi.insert(voronoi.end(), temp.begin(), temp.end());
    voronoi.insert(voronoi.end(), mid_edge.begin(), mid_edge.end());

    for(int i = 0; i < voronoi.size(); i++) {
        bool b = 0;
        for(int j = 0; j < voronoi.size(); j++) {
            if(i == j) continue;
            if(voronoi[i].start == voronoi[j].start || voronoi[i].start == voronoi[j].end ||
               voronoi[i].end == voronoi[j].start || voronoi[i].end == voronoi[j].end)
                b = 1;
        }
        if(b == 0) voronoi.erase(voronoi.begin() + i--);
    }

    return voronoi;
}

void read_points(const string& filename) {
    ifstream infile(filename);
    double x, y;
    while (infile >> x >> y) {
        Point point = {x, y};
        if (find(points.begin(), points.end(), point) == points.end()) {
            points.push_back(point);
        }
    }
}

void write_edges(const string& filename) {
    ofstream outfile(filename);
    for (const auto& edge : edges) {
        outfile << edge.start.x << " " << edge.start.y << " " << edge.end.x << " " << edge.end.y << "\n";
    }
    outfile << "\n";
    outfile.close();
    outfile.open("hyperplane/hyperplane" + to_string(file_num) + ".txt");
    outfile.close();
    outfile.open("convexhull/convexhull" + to_string(file_num++) + ".txt");
    for(auto i : convexhull_point) {
        outfile << "L " << i.x << ' ' << i.y << endl;
    }
    outfile.close();
}

int main() {
    read_points("points.txt");
    sort(points.begin(), points.end());
    edges = recursiveVoronoi(0, points.size());
    write_edges("lines.txt");
    return 0;
}
